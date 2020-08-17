#include "Output_writer.h"

#include <fmt/format.h>
#include <H5Cpp.h>
#include "array_utils.h"

void write_array(const array3d & array, const std::string name, H5::H5File file) {
    const hsize_t n1 = static_cast<hsize_t>(array.get_n1());
    const hsize_t n2 = static_cast<hsize_t>(array.get_n2());
    const hsize_t n3 = static_cast<hsize_t>(array.get_n3());
    const hsize_t dims[3] {n1, n2, n3};
    H5::DataSpace dataspace(3, dims);

    H5::DataSet dataset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(array(0,0,0)), H5::PredType::NATIVE_DOUBLE);

    const hsize_t attribute_dims[1] {3};
    H5::DataSpace attribute_data_space(1, attribute_dims);
    auto attribute = dataset.createAttribute("steps", H5::PredType::NATIVE_DOUBLE, attribute_data_space);
    const auto steps = array.get_steps();
    attribute.write(H5::PredType::NATIVE_DOUBLE, &steps);

    attribute = dataset.createAttribute("origin", H5::PredType::NATIVE_DOUBLE, attribute_data_space);
    const auto origin = array.get_origin();
    attribute.write(H5::PredType::NATIVE_DOUBLE, &origin);
}

void write_array(const array2d & array, const std::string name, H5::H5File file) {
    const hsize_t n1 = static_cast<hsize_t>(array.get_n1());
    const hsize_t n2 = static_cast<hsize_t>(array.get_n2());
    
    const hsize_t dims[2] {n1, n2};
    H5::DataSpace dataspace(2, dims);

    H5::DataSet dataset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(array(0,0)), H5::PredType::NATIVE_DOUBLE);

    const hsize_t attribute_dims[1] {2};
    H5::DataSpace attribute_data_space(1, attribute_dims);
    auto attribute = dataset.createAttribute("steps", H5::PredType::NATIVE_DOUBLE, attribute_data_space);
    const auto steps = array.get_steps();
    attribute.write(H5::PredType::NATIVE_DOUBLE, &(steps[0]));

    attribute = dataset.createAttribute("origin", H5::PredType::NATIVE_DOUBLE, attribute_data_space);
    const auto origin = array.get_origin_2d();
    attribute.write(H5::PredType::NATIVE_DOUBLE, &(origin[0]));

    auto plane = array.get_plane();
    if (plane != Plane::NONE) {
        H5::DataSpace scalar_data_space{};

        attribute = dataset.createAttribute("plane", H5::PredType::NATIVE_INT, scalar_data_space);
        auto plane_id = static_cast<int>(plane);
        attribute.write(H5::PredType::NATIVE_INT, &plane_id);

        attribute = dataset.createAttribute("plane_coordinate", H5::PredType::NATIVE_DOUBLE, scalar_data_space);
        auto plane_coordinate = array.get_plane_coordinate();
        attribute.write(H5::PredType::NATIVE_DOUBLE, &plane_coordinate);
    }
}

void initialize_slice_array(ivector3d size, vector3d steps, vector3d origin, const std::string name, H5::H5File file) {
    const hsize_t n1 = static_cast<hsize_t>(size.x);
    const hsize_t n2 = static_cast<hsize_t>(size.y);
    const hsize_t n3 = static_cast<hsize_t>(size.z);
    
    const hsize_t dims[3] {n1, n2, n3};
    H5::DataSpace dataspace(3, dims);

    H5::DataSet dataset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);

    const hsize_t attribute_dims[1] {3};
    H5::DataSpace attribute_data_space(1, attribute_dims);
    auto attribute = dataset.createAttribute("steps", H5::PredType::NATIVE_DOUBLE, attribute_data_space);
    attribute.write(H5::PredType::NATIVE_DOUBLE, &steps);

    attribute = dataset.createAttribute("origin", H5::PredType::NATIVE_DOUBLE, attribute_data_space);
    attribute.write(H5::PredType::NATIVE_DOUBLE, &origin);
}

void initialize_slice_array(ivector2d size, vector2d steps, vector2d origin, Plane plane, double plane_coordinate,
                            const std::string name, H5::H5File file) {
    const hsize_t n1 = static_cast<hsize_t>(size[0]);
    const hsize_t n2 = static_cast<hsize_t>(size[1]);
    
    const hsize_t dims[2] {n1, n2};
    H5::DataSpace dataspace(2, dims);

    H5::DataSet dataset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);

    const hsize_t attribute_dims[1] {2};
    H5::DataSpace attribute_data_space(1, attribute_dims);
    auto attribute = dataset.createAttribute("steps", H5::PredType::NATIVE_DOUBLE, attribute_data_space);
    attribute.write(H5::PredType::NATIVE_DOUBLE, &(steps[0]));

    attribute = dataset.createAttribute("origin", H5::PredType::NATIVE_DOUBLE, attribute_data_space);
    attribute.write(H5::PredType::NATIVE_DOUBLE, &(origin[0]));

    if (plane != Plane::NONE) {
        H5::DataSpace scalar_data_space{};

        attribute = dataset.createAttribute("plane", H5::PredType::NATIVE_INT, scalar_data_space);
        auto plane_id = static_cast<int>(plane);
        attribute.write(H5::PredType::NATIVE_INT, &plane_id);

        attribute = dataset.createAttribute("plane_coordinate", H5::PredType::NATIVE_DOUBLE, scalar_data_space);
        attribute.write(H5::PredType::NATIVE_DOUBLE, &plane_coordinate);
    }
}

void write_slice(const array2d & slice, int slice_index, const std::string name, H5::H5File file) {
    const hsize_t n1 = static_cast<hsize_t>(slice.get_n1());
    const hsize_t n2 = static_cast<hsize_t>(slice.get_n2());

    const hsize_t dims[2] {n1, n2};
    H5::DataSpace memory_dataspace(2, dims);
    
    auto dataset = file.openDataSet(name);
    auto write_dataspace = dataset.getSpace();
    const hsize_t count[3] {1, n1, n2};
    const hsize_t start[3] {static_cast<hsize_t>(slice_index), 0, 0};
    write_dataspace.selectHyperslab(H5S_SELECT_SET, count, start);

    dataset.write(&(slice(0, 0)), H5::PredType::NATIVE_DOUBLE, memory_dataspace, write_dataspace);
}

void write_slice(const array1d & slice, int slice_index, const std::string name, H5::H5File file) {
    const hsize_t n = static_cast<hsize_t>(slice.get_size());

    const hsize_t dims[1] {n};
    H5::DataSpace memory_dataspace(1, dims);
    
    auto dataset = file.openDataSet(name);
    auto write_dataspace = dataset.getSpace();
    const hsize_t count[2] {1, n};
    const hsize_t start[2] {static_cast<hsize_t>(slice_index), 0};
    write_dataspace.selectHyperslab(H5S_SELECT_SET, count, start);

    dataset.write(&(slice(0)), H5::PredType::NATIVE_DOUBLE, memory_dataspace, write_dataspace);
}

Output_writer::Output_writer(Output_parameters output_parameters, int count) :
    output_parameters(output_parameters) {
    if (output_parameters.output3d) {
        fields_file = H5::H5File(fmt::format("Fields_{:03d}.h5", count), H5F_ACC_TRUNC);
    }
    if (output_parameters.output_xy) {
        fields_xy_file = H5::H5File(fmt::format("Fields_xy_{:03d}.h5", count), H5F_ACC_TRUNC);
    }
    if (output_parameters.output_bunch) {
        bunch_parameters_file = H5::H5File(fmt::format("Bunch_{:03d}.h5", count), H5F_ACC_TRUNC);
    }
}

void Output_writer::initialize_slice_array(ivector3d size, vector3d steps, const array2d & array,
                                           const std::string name) {
    if (output_parameters.output3d) {
        ::initialize_slice_array(size, steps, array.get_origin_3d(), name, fields_file);
    }
    if (output_parameters.output_xy) {
        const ivector2d size_2d = {size.x, size.y};
        const vector2d steps_2d = {steps.x, steps.y};
        const vector2d origin_2d = {array.get_origin_3d().x, array.get_origin_3d().y};
        
        ::initialize_slice_array(size_2d, steps_2d, origin_2d, Plane::XY, output_parameters.z0, name, fields_xy_file);
    }
}

void Output_writer::write_array(array3d & array, std::string name) {
    if (output_parameters.output3d) {
        ::write_array(array, name, fields_file);
    }
    if (output_parameters.output_xy) {
        ::write_array(calculate_xy_slice(array, output_parameters.z0), name, fields_xy_file);
    }
}

void Output_writer::write_slice(array2d & array, std::string name, int slice_index) {
    if (output_parameters.output3d) {
        ::write_slice(array, slice_index, name, fields_file);
    }
    if (output_parameters.output_xy) {
        ::write_slice(calculate_y_slice(array, output_parameters.z0), slice_index, name, fields_xy_file);
    }
}

void Output_writer::write_bunch_parameters(const std::vector<bunch_particle_3d> & particles) {
    if (!output_parameters.output_bunch) {
        return;
    }
    
    const hsize_t size = particles.size();
    const hsize_t dims[1] {size};

    const hsize_t scalar_dims[1] {1};
    H5::DataSpace scalar_data_space{1, scalar_dims};

    H5::DataSpace dataspace(1, dims);


    H5::DataSet x_dataset = bunch_parameters_file.createDataSet("x", H5::PredType::NATIVE_DOUBLE, dataspace);
    H5::DataSet y_dataset = bunch_parameters_file.createDataSet("y", H5::PredType::NATIVE_DOUBLE, dataspace);
    H5::DataSet z_dataset = bunch_parameters_file.createDataSet("z", H5::PredType::NATIVE_DOUBLE, dataspace);
    H5::DataSet px_dataset = bunch_parameters_file.createDataSet("px", H5::PredType::NATIVE_DOUBLE, dataspace);
    H5::DataSet py_dataset = bunch_parameters_file.createDataSet("py", H5::PredType::NATIVE_DOUBLE, dataspace);
    H5::DataSet pz_dataset = bunch_parameters_file.createDataSet("pz", H5::PredType::NATIVE_DOUBLE, dataspace);
    H5::DataSet gamma_dataset = bunch_parameters_file.createDataSet("gamma", H5::PredType::NATIVE_DOUBLE, dataspace);
    H5::DataSet ex_dataset = bunch_parameters_file.createDataSet("ex", H5::PredType::NATIVE_DOUBLE, dataspace);
    H5::DataSet ey_dataset = bunch_parameters_file.createDataSet("ey", H5::PredType::NATIVE_DOUBLE, dataspace);
    H5::DataSet ez_dataset = bunch_parameters_file.createDataSet("ez", H5::PredType::NATIVE_DOUBLE, dataspace);
    H5::DataSet by_dataset = bunch_parameters_file.createDataSet("by", H5::PredType::NATIVE_DOUBLE, dataspace);
    H5::DataSet bz_dataset = bunch_parameters_file.createDataSet("bz", H5::PredType::NATIVE_DOUBLE, dataspace);

    for (hsize_t i = 0; i < size; i++) {
        hsize_t coords[1] = {i};

        dataspace.selectElements(H5S_SELECT_SET, 1, coords);
        const auto & p = particles[i];
        x_dataset.write(&(p.x), H5::PredType::NATIVE_DOUBLE, scalar_data_space, dataspace);
        y_dataset.write(&(p.y), H5::PredType::NATIVE_DOUBLE, scalar_data_space, dataspace);
        z_dataset.write(&(p.z), H5::PredType::NATIVE_DOUBLE, scalar_data_space, dataspace);
        px_dataset.write(&(p.px), H5::PredType::NATIVE_DOUBLE, scalar_data_space, dataspace);
        py_dataset.write(&(p.py), H5::PredType::NATIVE_DOUBLE, scalar_data_space, dataspace);
        pz_dataset.write(&(p.pz), H5::PredType::NATIVE_DOUBLE, scalar_data_space, dataspace);
        gamma_dataset.write(&(p.gamma), H5::PredType::NATIVE_DOUBLE, scalar_data_space, dataspace);
        ex_dataset.write(&(p.ex), H5::PredType::NATIVE_DOUBLE, scalar_data_space, dataspace);
        ey_dataset.write(&(p.ey), H5::PredType::NATIVE_DOUBLE, scalar_data_space, dataspace);
        ez_dataset.write(&(p.ez), H5::PredType::NATIVE_DOUBLE, scalar_data_space, dataspace);
        by_dataset.write(&(p.by), H5::PredType::NATIVE_DOUBLE, scalar_data_space, dataspace);
        bz_dataset.write(&(p.bz), H5::PredType::NATIVE_DOUBLE, scalar_data_space, dataspace);
    }
}