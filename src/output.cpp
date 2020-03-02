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
}