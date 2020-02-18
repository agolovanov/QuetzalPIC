#include "output.h"

void write_array(const array3d & array, const std::string name, H5::H5File file) {
    const hsize_t dims[3] {array.get_n1(), array.get_n2(), array.get_n3()};
    H5::DataSpace dataspace(3, dims);

    H5::DataSet dataset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(array(0,0,0)), H5::PredType::NATIVE_DOUBLE);
}