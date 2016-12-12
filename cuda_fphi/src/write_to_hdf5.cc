
#include <H5Cpp.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <string>
#include <iostream>
using namespace H5;

int write_data_to_hdf5(const char * out_filename,  std::string labels [] , float * data, size_t dim){
	try{
		size_t dst_size =  sizeof( float);
		size_t dst_offset[dim];
		hsize_t dims[1] = {dim};

		hid_t      field_type[dim];
		float fill_data[dim];







	 	 H5File file_out(out_filename, H5F_ACC_TRUNC);




	 	 DataSpace dataspace(1, dims);

	 	 hsize_t label_chunk_dims[1] = {dim};
	 	 DSetCreatPropList  plist_label;
	 	 plist_label.setChunk(1, label_chunk_dims);

	 	 plist_label.setDeflate(9);

	 	 H5::StrType datatype(H5::PredType::C_S1, H5T_VARIABLE);

	 	 DataSet data_out = file_out.createDataSet("/headers", datatype, dataspace, plist_label);

	 	 data_out.write(labels, datatype, dataspace);
	 	 plist_label.close();

	 	 dataspace.close();
	 	 data_out.close();


	 	 hsize_t data_dims[2] = {dim, dim};
	 	 hsize_t chunk_dims[2] = {1, dim};

	 	 DataSpace data_space(2, data_dims);
	 	 DSetCreatPropList  plist;
	 	 plist.setChunk(2, chunk_dims);

	 	 plist.setDeflate(9);


	 	 DataSet data_out_2 = file_out.createDataSet("/matrix", PredType::NATIVE_FLOAT, data_space, plist);

	 	 data_out_2.write(data, PredType::NATIVE_FLOAT, dataspace);


	 	 data_out_2.close();
	 	 plist.close();

	 	 data_space.close();

	 	 file_out.close();
	}
    catch(FileIException error)
    {
    	std::cout << "Error reading or writing hdf5 file " << std::string(out_filename) << "\n";
    	error.printError();
		return 0;
    }

    catch(DataSetIException error)
    {
    	std::cout << "Error reading or writing hdf5 file " << std::string(out_filename) << "\n";

    	error.printError();
		return 0;
    }

    catch(DataSpaceIException error)
    {
    	std::cout << "Error reading or writing hdf5 file " << std::string(out_filename) << "\n";

    	error.printError();
		return 0;
    }


    catch(...)
    {
    	std::cout << "Unknown error writing hdf5 file " << std::string(out_filename) << "\n";

    	return 0;
    }


    return 1;


}


int write_to_hdf5_qsub(const char * out_filename, float * data, size_t qsub_batch_size, size_t qsub_shared_batch_size, size_t node_number, size_t dim){
	try{
		H5File file_in (out_filename, H5F_ACC_RDWR);

		DataSet dataset = file_in.openDataSet("/matrix");

		DataSpace dataspace = dataset.getSpace();

		hsize_t dimsext[2] = {1, dim};
		hsize_t offset[2] = {qsub_shared_batch_size*node_number, 0};

		if(qsub_shared_batch_size == dim){
			dimsext[0] = 1;
			dimsext[1] =  dim;
			for(hsize_t row = 0 ; row < dim ; row++){
				offset[0] = row;
				offset[1] = 0;
				dataspace.selectHyperslab(H5S_SELECT_SET, dimsext, offset);

				DataSpace memspace(2, dimsext, NULL);

				dataset.write(&data[row*dim], PredType::NATIVE_FLOAT, memspace, dataspace);

				memspace.close();

			}

			dataspace.close();

			dataset.close();

			file_in.close();
		}else{


			for(hsize_t row = 0 ; row < qsub_batch_size; row++){
				offset[0] = qsub_shared_batch_size*node_number + row;
				offset[1] = 0;
				DataSpace memspace(2, dimsext, NULL);

				dataspace.selectHyperslab(H5S_SELECT_SET, dimsext, offset);


				dataset.write(&data[row*dim], PredType::NATIVE_FLOAT, memspace, dataspace);

				memspace.close();


			}



			dataspace.close();

			dataset.close();

			file_in.close();

		}

	}



    catch(FileIException error)
    {
    	std::cout << "Error reading or writing hdf5 file " << std::string(out_filename) << "\n";
    	error.printError();
		return 0;
    }

    catch(DataSetIException error)
    {
    	std::cout << "Error reading or writing hdf5 file " << std::string(out_filename) << "\n";

    	error.printError();
		return 0;
    }

    catch(DataSpaceIException error)
    {
    	std::cout << "Error reading or writing hdf5 file " << std::string(out_filename) << "\n";

    	error.printError();
		return 0;
    }


    catch(...)
    {
    	std::cout << "Unknown error writing hdf5 file " << std::string(out_filename) << "\n";

    	return 0;
    }


    return 1;



}



