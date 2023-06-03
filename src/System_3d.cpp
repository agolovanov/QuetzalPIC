#include <cmath>
#include <iostream>
#include <stdexcept>
#include <H5Cpp.h>
#include <cassert>
#include <fmt/format.h>
#include <vector>
#include <string>

#include "System_3d.h"
#include "containers.h"
#include "array2d.h"
#include "array3d.h"
#include "array_utils.h"
#include "Output_writer.h"
#include "constants.h"

std::string memory_formatter(long bytes) {
    double res = bytes;
    const std::vector<std::string> prefixes {"  B", "KiB", "MiB", "GiB", "TiB", "PiB"};
    int index = 0;
    while (res > 1024) {
        res /= 1024;
        index++;
    }
    return fmt::format("{:5.1f} {}", res, prefixes[index]);
}

void normalize_potential(array2d & psi) {
    int n1 = psi.get_n1();
    int n2 = psi.get_n2();
    double border_value = 0.0;
    for (int i = 0; i < n1; i++) {
        border_value += psi(i, 0);
    }
    for (int j = 1; j < n2; j++) {
        border_value += psi(0, j);
    }
    border_value /= (n1 + n2 - 1);
    #pragma omp parallel for
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            psi(i, j) -= border_value;
        }
    }
}

const std::string ASQR = "aSqr";
const std::string PSI = "psi";
const std::string EX = "ex";
const std::string EY = "ey";
const std::string EZ = "ez";
const std::string BY = "by";
const std::string BZ = "bz";
const std::string RHO = "rho";
const std::string JX = "jx";
const std::string JY = "jy";
const std::string JZ = "jz";
const std::string SUSCEPTIBILITY = "susceptibility";
const std::string RHO_BUNCH = "rho_bunch";
const std::string JX_BUNCH = "jx_bunch";
const std::string PARTICLE_ENERGY_DENSITY = "w_particle";
const std::string PARTICLE_SX = "sx_particle";
const std::string EM_ENERGY_DENSITY = "w_em";
const std::string EM_SX = "sx_em";

System_3d::System_3d(System_parameters & params, std::ostream & out) : 
    l(params.l),
    dt(params.dt),
    ppcy(params.ppcy),
    ppcz(params.ppcz),
    plasma_profile(params.plasma_profile),
    magnetic_field_iterations(params.magnetic_field_iterations),
    plasma_units(params.base_frequency_SI),
    output_parameters(params.output_parameters),
    out(out),
    species(params.species)
{
    if (magnetic_field_iterations < 0) {
        throw std::invalid_argument("Magnetic field iterations should be non-negative");
    }
    if (params.ppcy <= 0) {
        throw std::invalid_argument("ppcy should be positive");
    }
    if (params.ppcz <= 0) {
        throw std::invalid_argument("ppcz should be positive");
    }

    n.x = static_cast<int>(l.x / params.d.x);
    n.y = static_cast<int>(l.y / params.d.y);
    n.z = static_cast<int>(l.z / params.d.z);

    out << fmt::format("Spatial dimensions:  [{}, {}, {}]", n.x, n.y, n.z) << std::endl;

    d.x = l.x / n.x;
    d.y = l.y / n.y;
    d.z = l.z / n.z;

    out << fmt::format("Steps:               [{}, {}, {}]", d.x, d.y, d.z) << std::endl;
    out << fmt::format("Simulation box size: [{}, {}, {}]", l.x, l.y, l.z) << std::endl;

    time_iterations = static_cast<int>(params.t_end / dt) + 1;

    t_end = dt * (time_iterations - 1);

    out << fmt::format("Timestep: {}, end time: {}, iterations: {}", dt, t_end, time_iterations) << std::endl;

    out << "----------------------------------------" << std::endl;
    out << "PLASMA UNITS\n";
    out << fmt::format("Frequency {:.4g} rad/s\n", plasma_units.frequency);
    out << fmt::format("Density {:.4g} cm^-3\n", plasma_units.density);
    out << fmt::format("Wavelength {:.4g} cm\n", plasma_units.wavelength);
    out << fmt::format("Schwinger field in plasma units: {:.4g}\n", plasma_units.field_schwinger);
    
    out << "----------------------------------------" << std::endl;

    out << "REGISTERED SPECIES\n";
    for (auto & species : species.get_species()) {
        out << fmt::format("{}: charge {:.4g}, mass {:.4g}, charge_to_mass_ratio {:.4g}\n", species.name, species.charge, species.mass, species.charge_to_mass_ratio);
    }    
    out << "----------------------------------------" << std::endl;

    const int bunches_count = params.bunch_parameters_array.size();
    std::vector<size_t> bunch_particles_count_array(bunches_count);
    size_t bunch_particles_count = 0;
    std::vector<long> bunch_particle_memory_array(bunches_count);
    long bunch_particle_memory = 0;
    
    for (int i = 0; i < bunches_count; i++) {
        const auto & bunch = params.bunch_parameters_array[i];
        bunch_particles_count_array[i] = count_bunch_particles(bunch.ppc, bunch.rho);
        bunch_particles_count += bunch_particles_count_array[i];
        bunch_particle_memory_array[i] = static_cast<long>(sizeof(bunch_particle_3d)) * bunch_particles_count_array[i];
        bunch_particle_memory += bunch_particle_memory_array[i];
    }

    const size_t wake_particles_count = count_wake_particles(ppcy, ppcz, plasma_profile);

    const long fourier_memory = sizeof(double) * n.y * n.z + 2 * sizeof(double) * n.y * (n.z / 2 + 1);
    const long array2d_memory = 14l * sizeof(double) * n.y * n.z;
    const long array3d_memory = 11l * sizeof(double) * n.x * n.y * n.z;
    const long wake_particle_memory = static_cast<long>(sizeof(wake_particle_2d)) * wake_particles_count;
    const long total_memory = array2d_memory + array3d_memory + wake_particle_memory + bunch_particle_memory + fourier_memory;

    out << "Expected RAM usage:\n";
    out << fmt::format("3D arrays:       {}\n", memory_formatter(array3d_memory))
        << fmt::format("2D arrays:       {}\n", memory_formatter(array2d_memory))
        << fmt::format("Fourier:         {}\n", memory_formatter(fourier_memory))
        << fmt::format("Wake particles:  {}\n", memory_formatter(wake_particle_memory))
        << fmt::format("Bunch particles: {}\n", memory_formatter(bunch_particle_memory))
        << fmt::format("Total:           {}", memory_formatter(total_memory)) << std::endl;

    out << "----------------------------------------" << std::endl;


    bunch_particles = std::vector<bunch_particle_3d>{bunch_particles_count};

    size_t index = 0;
    for (int i = 0; i < bunches_count; i++) {
        const auto & bunch = params.bunch_parameters_array[i];
        init_bunch_particles(index, bunch);
        index += bunch_particles_count_array[i];
    }

    wake_particles = std::vector<wake_particle_2d>{wake_particles_count};
    
    fourier = Fourier2d(n.y, n.z);

    const ivector2d size_yz{n.y, n.z};
    const vector2d d_yz{d.y, d.z};
    
    psi_middle = array2d(size_yz, d_yz);
    djy_dxi = array2d(size_yz, d_yz, {0.5 * d.y, 0});
    djz_dxi = array2d(size_yz, d_yz, {0, 0.5 * d.z});
    rho_ion = array2d(size_yz, d_yz);

    a_sqr = array3d(n, d);
    susceptibility = array3d(n, d);
    rho_bunch = array3d(n, d);
    jx_bunch = array3d(n, d);
    by = array3d(n, d, {0, 0, 0.5 * d.z});
    bz = array3d(n, d, {0, 0.5 * d.y, 0});
    ex = array3d(n, d, {-0.5 * d.x, 0, 0});
    ey = array3d(n, d, {0, 0.5 * d.y, 0});
    ez = array3d(n, d, {0, 0, 0.5 * d.z});

    em_energy_density = array3d(n, d, {0, 0, 0});
    em_sx = array3d(n, d, {0, 0, 0});

    psi = array2d(size_yz, d, {0, 0, 0}, Plane::YZ);
    psi_prev = array2d(size_yz, d, {0, 0, 0}, Plane::YZ);
    rho = array2d(size_yz, d, {0, 0, 0}, Plane::YZ);
    jx = array2d(size_yz, d, {0, 0, 0}, Plane::YZ);
    jy = array2d(size_yz, d, {-0.5 * d.x, 0.5 * d.y, 0}, Plane::YZ);
    jz = array2d(size_yz, d, {-0.5 * d.x, 0, 0.5 * d.z}, Plane::YZ);
    jy_next = array2d(size_yz, d, {-0.5 * d.x, 0.5 * d.y, 0}, Plane::YZ);
    jz_next = array2d(size_yz, d, {-0.5 * d.x, 0, 0.5 * d.z}, Plane::YZ);

    particle_energy_density = array2d(size_yz, d, {0, 0, 0}, Plane::YZ);
    particle_sx = array2d(size_yz, d, {0, 0, 0}, Plane::YZ);

    init_a_sqr(params.a_sqr);
}

void System_3d::run() {
    for (int ti = 0; ti < time_iterations; ti++) {
        out << "Iteration " << ti << std::endl;

        out << "Depositing densities..." << std::endl;

        #pragma omp parallel for
        for (int i = 0; i < n.x; i++) {
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    rho_bunch(i, j, k) = 0.0;
                    jx_bunch(i, j, k) = 0.0;
                }
            }
        }

        const int bunch_particles_count = bunch_particles.size();

        #pragma omp parallel for
        for (int pi = 0; pi < bunch_particles_count; pi++) {
            auto & p = bunch_particles[pi];
            const double charge_density = species[p.species_id].charge * p.n;
            deposit(p.x, p.y, p.z, charge_density, rho_bunch);
            deposit(p.x, p.y, p.z, charge_density * p.px / p.gamma, jx_bunch);
        }

        out << "Solving wakefield..." << std::endl;

        solve_wakefield(ti);

        

        if (ti < time_iterations - 1) {
            out << "Updating particles..." << std::endl;

            #pragma omp parallel for
            for (int pi = 0; pi < bunch_particles_count; pi++) {
                auto & p = bunch_particles[pi];
                const double cmr = species[p.species_id].charge_to_mass_ratio;
                
                const double fx = cmr * (p.ex + (p.py * p.bz - p.pz * p.by) / p.gamma);
                const double fy = cmr * (p.ey - p.px * p.bz / p.gamma);
                const double fz = cmr * (p.ez + p.px * p.by / p.gamma);

                p.x += dt * (p.px / p.gamma - 1);
                p.y += dt * p.py / p.gamma;
                p.z += dt * p.pz / p.gamma;
                p.px += dt * fx;
                p.py += dt * fy;
                p.pz += dt * fz;
                p.gamma = sqrt(1 + p.px * p.px + p.py * p.py + p.pz * p.pz);
            }
        }
    }
}

void System_3d::solve_wakefield(int iteration) {
    auto output_writer = Output_writer(output_parameters, iteration);

    std::vector<Output_reference<array3d>> output_arrays_3d;
    output_arrays_3d.push_back(Output_reference<array3d>(ASQR, &a_sqr));
    output_arrays_3d.push_back(Output_reference<array3d>(SUSCEPTIBILITY, &susceptibility));
    output_arrays_3d.push_back(Output_reference<array3d>(RHO_BUNCH, &rho_bunch));
    output_arrays_3d.push_back(Output_reference<array3d>(JX_BUNCH, &jx_bunch));
    output_arrays_3d.push_back(Output_reference<array3d>(EX, &ex));
    output_arrays_3d.push_back(Output_reference<array3d>(EY, &ey));
    output_arrays_3d.push_back(Output_reference<array3d>(EZ, &ez));
    output_arrays_3d.push_back(Output_reference<array3d>(BY, &by));
    output_arrays_3d.push_back(Output_reference<array3d>(BZ, &bz));
    output_arrays_3d.push_back(Output_reference<array3d>(EM_ENERGY_DENSITY, &em_energy_density));
    output_arrays_3d.push_back(Output_reference<array3d>(EM_SX, &em_sx));

    std::vector<Output_reference<array2d>> output_arrays_2d;
    output_arrays_2d.push_back(Output_reference<array2d>(PSI, &psi));
    output_arrays_2d.push_back(Output_reference<array2d>(RHO, &rho));
    output_arrays_2d.push_back(Output_reference<array2d>(JX, &jx));
    output_arrays_2d.push_back(Output_reference<array2d>(JY, &jy));
    output_arrays_2d.push_back(Output_reference<array2d>(JZ, &jz));
    output_arrays_2d.push_back(Output_reference<array2d>(PARTICLE_ENERGY_DENSITY, &particle_energy_density));
    output_arrays_2d.push_back(Output_reference<array2d>(PARTICLE_SX, &particle_sx));

    for (auto & output_arr : output_arrays_2d) {
        output_writer.initialize_slice_array(n, d, *(output_arr.ptr), output_arr.name);
    }

    init_wake_particles(ppcy, ppcz, plasma_profile);
    int particle_number = wake_particles.size();

    #pragma omp parallel for
    for (int j = 0; j < n.y; j++) {
        for (int k = 0; k < n.z; k++) {
            ex(n.x - 1, j, k) = 0;
        }
    }
    
    #pragma omp parallel for
    for (int j = 0; j < n.y; j++) {
        for (int k = 0; k < n.z; k++) {
            psi(j, k) = 0;
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < n.x; i++) {
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                by(i, j, k) = 0;
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < n.x; i++) {
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                bz(i, j, k) = 0;
            }
        }
    }

    #pragma omp parallel for
    for (int j = 0; j < n.y; j++) {
        for (int k = 0; k < n.z; k++) {
            rho_ion(j, k) = 0;
        }
    }

    #pragma omp parallel for
    for (int j = 0; j < n.y; j++) {
        for (int k = 0; k < n.z; k++) {
            jy(j, k) = 0;
        }
    }

    #pragma omp parallel for
    for (int j = 0; j < n.y; j++) {
        for (int k = 0; k < n.z; k++) {
            jz(j, k) = 0;
        }
    }

    // rho_ion deposition
    #pragma omp parallel for
    for (int pi = 0; pi < particle_number; pi++) {
        auto & p = wake_particles[pi];
        deposit(p.y, p.z, -p.n, rho_ion);
    }

    bool stop_flag = false;
    const double THRESHOLD_B = 100;

    for (int i = n.x - 1; i >= 0; i--) {

        out << "Slice " << i << std::endl;

        #pragma omp parallel for
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                psi_prev(j, k) = psi(j, k);
            }
        }

        // psi source deposition
        #pragma omp parallel for
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                fourier.in[n.z * j + k] = jx_bunch(i, j, k) - rho_bunch(i, j, k) - rho_ion(j, k);
            }
        }

        #pragma omp parallel for
        for (int pi = 0; pi < particle_number; pi++) {
            auto & p = wake_particles[pi];
            deposit(p.y, p.z, -p.n, fourier.in, {d.y, d.z}, {n.y, n.z});
        }

        // cacluate psi
        solve_poisson_equation();

        #pragma omp parallel for
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                psi(j, k) = (fourier.in[n.z * j + k]) / n.y / n.z;
            }
        }
        normalize_potential(psi);
        increase_minimum(psi, psi_threshold - 1);

        // calculate initial gamma, px
        #pragma omp parallel for
        for (int pi = 0; pi < particle_number; pi++) {
            auto & p = wake_particles[pi];
            double a_particle = array_to_particle(p.y, p.z, a_sqr, i);
            double psi_particle = array_to_particle(p.y, p.z, psi);
            p.gamma = 0.5 * (1 + p.py * p.py + p.pz * p.pz + a_particle + (1 + psi_particle) * (1 + psi_particle))
                    / (1 + psi_particle);
            p.px = 0.5 * (1 + p.py * p.py + p.pz * p.pz + a_particle - (1 + psi_particle) * (1 + psi_particle))
                    / (1 + psi_particle);
            assert(p.gamma >= 1.0);
        }

        // initial jx, rho deposition
        #pragma omp parallel for
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                rho(j, k) = rho_bunch(i, j, k);
                jx(j, k) = jx_bunch(i, j, k);
                susceptibility(i, j, k) = 0.0;
            }
        }

        #pragma omp parallel for
        for (int pi = 0; pi < particle_number; pi++) {
            auto & p = wake_particles[pi];
            double vx = p.px / p.gamma;
            deposit(p.y, p.z, p.n * vx / (1 - vx), jx);
            deposit(p.y, p.z, p.n / (1 - vx), rho);
            deposit(p.y, p.z, -p.n / (1 - vx) / p.gamma, susceptibility, i);
        }

        if (i < n.x - 1) {
            #pragma omp parallel for
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    by(i, j, k) = by(i+1, j, k);
                    bz(i, j, k) = bz(i+1, j, k);
                }
            }
        }

        if (i == 0) {
            break;
        }

        for (int iteration = 0; iteration < magnetic_field_iterations; iteration++) {

            // advance momenta
            #pragma omp parallel for
            for (int pi = 0; pi < particle_number; pi++) {
                auto & p = wake_particles[pi];
                double psi_particle = array_to_particle(p.y, p.z, psi);
                double da_dy_particle = array_yder_to_particle(p.y, p.z, a_sqr, i);
                double da_dz_particle = array_zder_to_particle(p.y, p.z, a_sqr, i);
                double dpsi_dy_particle = array_yder_to_particle(p.y, p.z, psi);
                double dpsi_dz_particle = array_zder_to_particle(p.y, p.z, psi);
                double by_particle = array_to_particle(p.y, p.z, by, i);
                double bz_particle = array_to_particle(p.y, p.z, bz, i);

                p.py_next = p.py - d.x * 0.5 * da_dy_particle / (1 + psi_particle);
                p.py_next += d.x * p.gamma * dpsi_dy_particle / (1 + psi_particle);
                p.py_next -= d.x * bz_particle;
                p.pz_next = p.pz - d.x * 0.5 * da_dz_particle / (1 + psi_particle);
                p.pz_next += d.x * p.gamma * dpsi_dz_particle / (1 + psi_particle);
                p.pz_next += d.x * by_particle;
            }

            // advance half coordinate
            #pragma omp parallel for
            for (int pi = 0; pi < particle_number; pi++) {
                auto & p = wake_particles[pi];
                double psi_particle = array_to_particle(p.y, p.z, psi);
                p.y_middle = p.y + 0.5 * d.x * p.py_next / (1 + psi_particle);
                p.z_middle = p.z + 0.5 * d.x * p.pz_next / (1 + psi_particle);
                normalize_coordinates(p.y_middle, p.z_middle);
            }

            // psi_source middle deposition
            #pragma omp parallel for
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    fourier.in[n.z * j + k] = jx_bunch(i, j, k) - rho_bunch(i, j, k) - rho_ion(j, k);
                }
            }

            #pragma omp parallel for
            for (int pi = 0; pi < particle_number; pi++) {
                auto & p = wake_particles[pi];
                deposit(p.y_middle, p.z_middle, -p.n, fourier.in, {d.y, d.z}, {n.y, n.z});
            }

            // calculate psi_middle
            solve_poisson_equation();

            #pragma omp parallel for
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    psi_middle(j, k) = (fourier.in[n.z * j + k]) / n.y / n.z;
                }
            }
            normalize_potential(psi_middle);
            increase_minimum(psi_middle, psi_threshold - 1);

            // deposit jy, jz
            #pragma omp parallel for
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++){
                    jy_next(j, k) = 0.0;
                    jz_next(j, k) = 0.0;
                }
            }

            #pragma omp parallel for
            for (int pi = 0; pi < particle_number; pi++) {
                auto & p = wake_particles[pi];
                double psi_particle = array_to_particle(p.y_middle, p.z_middle, psi_middle);
                deposit(p.y_middle, p.z_middle, p.n * p.py_next / (1 + psi_particle), jy_next);
                deposit(p.y_middle, p.z_middle, p.n * p.pz_next / (1 + psi_particle), jz_next);
            }

            // calculate new djy_dxi, djz_dxi
            #pragma omp parallel for
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    djy_dxi(j, k) = (jy(j, k) - jy_next(j, k)) / d.x;
                    djz_dxi(j, k) = (jz(j, k) - jz_next(j, k)) / d.x;
                }
            }

            // calculate new gamma, px
            #pragma omp parallel for
            for (int pi = 0; pi < particle_number; pi++) {
                auto & p = wake_particles[pi];
                double a_particle = array_to_particle(p.y, p.z, a_sqr, i);
                double psi_particle = array_to_particle(p.y, p.z, psi);
                double p_squared = 0.25 * (p.py + p.py_next) * (p.py + p.py_next) + 0.25 * (p.pz + p.pz_next) * (p.pz + p.pz_next);
                p.gamma = 0.5 * (1 + p_squared + a_particle +
                        (1 + psi_particle) * (1 + psi_particle)) / (1 + psi_particle);
                p.px = 0.5 * (1 + p_squared + a_particle - (1 + psi_particle) * (1 + psi_particle)) / (1 + psi_particle);
                assert(p.gamma >= 1.0);
            }

            // new jx, rho deposition
            #pragma omp parallel for
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    rho(j, k) = rho_bunch(i, j, k);
                    jx(j, k) = jx_bunch(i, j, k);
                    susceptibility(i, j, k) = 0.0;
                }
            }

            #pragma omp parallel for
            for (int pi = 0; pi < particle_number; pi++) {
                auto & p = wake_particles[pi];
                double vx = p.px / p.gamma;
                deposit(p.y, p.z, p.n * vx / (1 - vx), jx);
                deposit(p.y, p.z, p.n / (1 - vx), rho);
                deposit(p.y, p.z, -p.n / (1 - vx) / p.gamma, susceptibility, i);
            }

            // new source for B_y
            #pragma omp parallel for
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z-1; k++) {
                    fourier.in[n.z * j + k] = - magnetic_field_D * by(i, j, k) - (jx(j, k+1) - jx(j, k)) / d.z + djz_dxi(j, k);
                }
                fourier.in[n.z * j + (n.z-1)] = - magnetic_field_D * by(i, j, n.z-1) - (jx(j, 0) - jx(j, n.z-1)) / d.z + djz_dxi(j, n.z-1);
            }

            // new guess for B_y
            solve_poisson_equation(magnetic_field_D);

            #pragma omp parallel for shared(stop_flag)
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    by(i, j, k) = fourier.in[n.z * j + k] / n.y / n.z;
                    if (fabs(by(i, j, k)) > THRESHOLD_B) {
                        stop_flag = true;
                    }
                    if (std::isnan(by(i, j, k))) {
                        stop_flag = true;
                    }
                }
            }

            // new source for B_z
            #pragma omp parallel for
            for (int j = 0; j < n.y-1; j++) {
                for (int k = 0; k < n.z; k++) {
                    fourier.in[n.z * j + k] = - magnetic_field_D * bz(i, j, k) + (jx(j+1, k) - jx(j, k)) / d.y - djy_dxi(j, k);
                }
            }
            #pragma omp parallel for
            for (int k = 0; k < n.z; k++) {
                fourier.in[n.z * (n.y-1) + k] = - magnetic_field_D * bz(i, 0, k) + (jx(0, k) - jx(n.y-1, k)) / d.y - djy_dxi(n.y-1, k);
            }

            // new guess for B_z
            solve_poisson_equation(magnetic_field_D);

            #pragma omp parallel for shared(stop_flag)
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    bz(i, j, k) = fourier.in[n.z * j + k] / n.y / n.z;
                    if (fabs(bz(i, j, k)) > THRESHOLD_B) {
                        stop_flag = true;
                    }
                    if (std::isnan(bz(i, j, k))) {
                        stop_flag = true;
                    }
                }
            }

            if (stop_flag) {
                break;
            }

        }

        if (stop_flag) {
            out << "Stop condition encountered, stopping the simulation" << std::endl;
            break;
        }

        #pragma omp parallel for
        for (int pi = 0; pi < particle_number; pi++) {
            auto & p = wake_particles[pi];
            p.py = p.py_next;
            p.pz = p.pz_next;
        }

        // advance coordinate
        #pragma omp parallel for
        for (int pi = 0; pi < particle_number; pi++) {
            auto & p = wake_particles[pi];
            double psi_particle = array_to_particle(p.y, p.z, psi_middle);
            p.y += d.x * p.py / (1 + psi_particle);
            p.z += d.x * p.pz / (1 + psi_particle);
            normalize_coordinates(p.y, p.z);
        }

        output_step(output_writer, output_arrays_2d, i);

        #pragma omp parallel for
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                jy(j, k) = jy_next(j, k);
                jz(j, k) = jz_next(j, k);
            }
        }
    }

    output_step(output_writer, output_arrays_2d, 0);

    calculate_em_density();

    for (auto & output_arr : output_arrays_3d) {
        output_writer.write_array(*(output_arr.ptr), output_arr.name);
    }


    const int particles_size = bunch_particles.size();
    #pragma omp parallel for
    for (int i = 0; i < particles_size; i++) {
        auto & p = bunch_particles[i];
        p.ex = array_to_particle(p.x, p.y, p.z, ex);
        p.ey = array_to_particle(p.x, p.y, p.z, ey);
        p.ez = array_to_particle(p.x, p.y, p.z, ez);
        p.by = array_to_particle(p.x, p.y, p.z, by);
        p.bz = array_to_particle(p.x, p.y, p.z, bz);

        const double fx = p.gamma * p.ex + p.py * p.bz - p.pz * p.by;
        const double fy = p.gamma * p.ey - p.px * p.bz;
        const double fz = p.gamma * p.ez + p.px * p.by;
        const double f_long = p.px * p.ex + p.py * p.ey + p.pz * p.ez;
        p.chi = sqrt(fx * fx + fy * fy + fz * fz - f_long * f_long) / plasma_units.field_schwinger;
    }

    output_writer.write_bunch_parameters(bunch_particles, d.x * d.y * d.z * plasma_units.number_density_norm);
}

void System_3d::output_step(Output_writer & output_writer, 
                            const std::vector<Output_reference<array2d>> & output_arrays_2d, int slice_index) {
    const int i = slice_index;
    
    // calculate ex
    if (i < n.x - 1) {
        #pragma omp parallel for
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                ex(i, j, k) = (psi(j, k) - psi_prev(j, k)) / d.x;
            }
        }
    }

    // calculate ey
    #pragma omp parallel for
    for (int j = 0; j < n.y-1; j++) {
        for (int k = 0; k < n.z; k++) {
            ey(i, j, k) = -(psi(j+1, k) - psi(j, k)) / d.y + bz(i, j, k);
        }
    }
    for (int k = 0; k < n.z; k++) {
        ey(i, n.y-1, k) = -(psi(0, k) - psi(n.y-1, k)) / d.y + bz(i, n.y-1, k);
    }
    
    // calculate ez
    #pragma omp parallel for
    for (int j = 0; j < n.y; j++) {
        for (int k = 0; k < n.z-1; k++) {
            ez(i, j, k) = -(psi(j, k+1) - psi(j, k)) / d.z - by(i, j, k);
        }
        ez(i, j, n.z-1) = -(psi(j, 0) - psi(j, n.z-1)) / d.z - by(i, j, n.z-1);
    }

    // deposit particle energy density
    #pragma omp parallel for
    for (int j = 0; j < n.y; j++) {
        for (int k = 0; k < n.z; k++) {
            particle_energy_density(j, k) = 0.0;
            particle_sx(j, k) = 0.0;
        }
    }

    int particle_number = wake_particles.size();

    #pragma omp parallel for
    for (int pi = 0; pi < particle_number; pi++) {
        auto & p = wake_particles[pi];
        double vx = p.px / p.gamma;
        deposit(p.y, p.z, - p.n * (p.gamma - 1) / (1 - vx), particle_energy_density);
        deposit(p.y, p.z, - p.n * (p.px - vx) / (1 - vx), particle_sx);
    }
    
    for (auto & output_arr : output_arrays_2d) {
        output_writer.write_slice(*(output_arr.ptr), output_arr.name, slice_index);
    }
}

void System_3d::normalize_coordinates(double & y, double & z) {
    y = fmod(y, l.y);
    if (y < 0) {
        y += l.y;
    }

    assert(y >= 0);
    assert(y < l.y);

    z = fmod(z, l.z);
    if (z < 0) {
        z += l.z;
    }
    
    assert(z >= 0);
    assert(z < l.z);
}

void System_3d::solve_poisson_equation(double D) {
    fourier.forward_transform();
    fourier.out[0][0] = 0;
    fourier.out[0][1] = 0;
    const int nz_size = (n.z / 2) + 1;

    #pragma omp parallel for
    for (int j = 0; j < n.y; j++) {
        for (int k = 0; k < nz_size; k++){
            if ((k == 0) && (j == 0)) {
                continue;
            }

            double multiplier = -D + (2 * cos(2 * M_PI * j / n.y) - 2) / d.y / d.y + (2 * cos(2 * M_PI * k / n.z) - 2) / d.z / d.z;

            fourier.out[nz_size * j + k][0] /= multiplier;
            fourier.out[nz_size * j + k][1] /= multiplier;
        }
    }
    fourier.backward_transform();
}

size_t System_3d::count_wake_particles(int ppcy, int ppcz, std::function<double(double, double)> plasma_profile) const {
    size_t count = 0;

    for (int i = 0; i < ppcy * n.y; i++) {
        for (int j = 0; j < ppcz * n.z; j++) {
            const double y = (i + 0.5) * d.y / ppcy;
            const double z = (j + 0.5) * d.z / ppcz;
            if (plasma_profile(y, z) != 0) {
                count++;
            }
        }
    }

    return count;
}

void System_3d::init_wake_particles(int ppcy, int ppcz, std::function<double(double, double)> plasma_profile) {
    assert(ppcy > 0);
    assert(ppcz > 0);

    size_t index = 0;

    for (int i = 0; i < ppcy * n.y; i++) {
        for (int j = 0; j < ppcz * n.z; j++) {
            const double y = (i + 0.5) * d.y / ppcy;
            const double z = (j + 0.5) * d.z / ppcz;
            double value = plasma_profile(y, z);
            if (value != 0) {
                wake_particles[index].y = y;
                wake_particles[index].z = z;
                wake_particles[index].py = 0;
                wake_particles[index].pz = 0;
                wake_particles[index].py_next = 0;
                wake_particles[index].pz_next = 0;
                wake_particles[index].gamma = 0;
                wake_particles[index].n = value / ppcy / ppcz;
                index++;
            }
        }
    }
}

size_t System_3d::count_bunch_particles(ivector3d ppc, std::function<double(double, double, double)> rho) const {
    assert(ppc.x > 0);
    assert(ppc.y > 0);
    assert(ppc.z > 0);
    
    size_t count = 0;
    for (int i = 0; i < ppc.x * n.x; i++) {
        for (int j = 0; j < ppc.y * n.y; j++) {
            for (int k = 0; k < ppc.z * n.z; k++) {
                const double x = (i + 0.5) * d.x / ppc.x;
                const double y = (j + 0.5) * d.y / ppc.y;
                const double z = (k + 0.5) * d.z / ppc.z;
                if (rho(x, y, z) > 0.0) {
                    count++;
                }
            }
        }
    }
    return count;
}

void System_3d::init_bunch_particles(size_t index, Bunch_parameters bunch) {
    assert(bunch.ppc.x > 0);
    assert(bunch.ppc.y > 0);
    assert(bunch.ppc.z > 0);
    
    for (int i = 0; i < bunch.ppc.x * n.x; i++) {
        for (int j = 0; j < bunch.ppc.y * n.y; j++) {
            for (int k = 0; k < bunch.ppc.z * n.z; k++) {
                const double x = (i + 0.5) * d.x / bunch.ppc.x;
                const double y = (j + 0.5) * d.y / bunch.ppc.y;
                const double z = (k + 0.5) * d.z / bunch.ppc.z;
                double value = bunch.rho(x, y, z);
                if (value > 0.0) {
                    bunch_particles[index].x = x;
                    bunch_particles[index].y = y;
                    bunch_particles[index].z = z;
                    bunch_particles[index].px = 0;
                    bunch_particles[index].py = 0;
                    bunch_particles[index].pz = 0;
                    bunch_particles[index].gamma = bunch.gamma;
                    bunch_particles[index].px = sqrt(bunch.gamma * bunch.gamma - 1);
                    bunch_particles[index].n = value / bunch.ppc.x / bunch.ppc.y / bunch.ppc.z;
                    bunch_particles[index].species_id = bunch.species_id;
                    index++;
                }
            }
        }
    }
}

void System_3d::init_a_sqr(std::function<double(double, double, double)> func) {
    for (int i = 0; i < n.x; i++) {
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                double x = i * d.x;
                double y = j * d.y;
                double z = k * d.z;
                a_sqr(i, j, k) = func(x, y, z);
            }
        }
    }
}

void System_3d::increase_minimum(array3d & array, int slice, double value) const {
    auto n = array.get_dimensions();
    #pragma omp parallel for
    for (int j = 0; j < n.y; j++) {
        for (int k = 0; k < n.z; k++) {
            if (array(slice, j, k) < value) {
                array(slice, j, k) = value;
            }
        }
    }
}

void System_3d::increase_minimum(array2d & array, double value) const {
    auto n1 = array.get_n1();
    auto n2 = array.get_n2();
    #pragma omp parallel for
    for (int j = 0; j < n1; j++) {
        for (int k = 0; k < n2; k++) {
            if (array(j, k) < value) {
                array(j, k) = value;
            }
        }
    }
}

void System_3d::calculate_em_density() {
    #pragma omp parallel for
    for (int i = 0; i < n.x; i++) {
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                double ex_point;
                if (i > 0) {
                    ex_point = 0.5 * (ex(i-1, j, k) + ex(i, j, k));
                } else {
                    ex_point = ex(i, j, k);
                }
                double ey_point = ey(i, j, k);
                double ez_point = ez(i, j, k);
                double by_point = by(i, j, k);
                double bz_point = bz(i, j, k);
                em_energy_density(i, j, k) = 0.5 * (ex_point * ex_point + ey_point * ey_point + ez_point * ez_point + by_point * by_point
                                                    + bz_point * bz_point);
                em_sx(i, j, k) = ey_point * bz_point - ez_point * by_point;
            }
        }
    }
}
