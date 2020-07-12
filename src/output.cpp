#include "output.h"

void write_array(const array3d & array, const std::string name, H5::H5File file) {
    const hsize_t dims[3] {array.get_n1(), array.get_n2(), array.get_n3()};
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
    const hsize_t dims[2] {array.get_n1(), array.get_n2()};
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
    const hsize_t dims[3] {size.x, size.y, size.z};
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
    const hsize_t dims[2] {size[0], size[1]};
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
    int n1 = slice.get_n1();
    int n2 = slice.get_n2();

    const hsize_t dims[2] {n1, n2};
    H5::DataSpace memory_dataspace(2, dims);
    
    auto dataset = file.openDataSet(name);
    auto write_dataspace = dataset.getSpace();
    const hsize_t count[3] {1, n1, n2};
    const hsize_t start[3] {slice_index, 0, 0};
    write_dataspace.selectHyperslab(H5S_SELECT_SET, count, start);

    dataset.write(&(slice(0, 0)), H5::PredType::NATIVE_DOUBLE, memory_dataspace, write_dataspace);
}

void write_slice(const array1d & slice, int slice_index, const std::string name, H5::H5File file) {
    int n = slice.get_size();

    const hsize_t dims[1] {n};
    H5::DataSpace memory_dataspace(1, dims);
    
    auto dataset = file.openDataSet(name);
    auto write_dataspace = dataset.getSpace();
    const hsize_t count[2] {1, n};
    const hsize_t start[2] {slice_index, 0};
    write_dataspace.selectHyperslab(H5S_SELECT_SET, count, start);

    dataset.write(&(slice(0)), H5::PredType::NATIVE_DOUBLE, memory_dataspace, write_dataspace);
}