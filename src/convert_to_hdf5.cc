#include "solar.h"
#include <hdf5.h>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include "safelib.h"
using namespace std;
static int convert_file(Tcl_Interp * interp, const char * filename, const char * output_filename){
	const char * errmsg = 0;
	SolarFile * file = SolarFile::open("hdf5 convert", filename,&errmsg); 

	if(errmsg){
		cout << errmsg << endl;
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
		

	int field_count;
	const char ** names = file->names(&field_count, &errmsg);
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}	
	file->start_setup(&errmsg);
	if(errmsg){
		RESULT_LIT(errmsg);
		free(names);
		return TCL_ERROR;
	}
	for(int field = 0 ; field < field_count; field++){
		file->setup(names[field], &errmsg);
		if(errmsg){
			RESULT_LIT(errmsg);
			free(names);
			return TCL_ERROR;
		}
	}		
	char ** file_data;
	string str;
	string dset_name; 
	hsize_t dims[1];
	hid_t file_id;
	file_id = H5Fcreate(output_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 
	hid_t datatype = H5Tcopy (H5T_C_S1);
	H5Tset_size (datatype, H5T_VARIABLE);
	dims[0] = field_count;		
	hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
	dset_name = "/header";
	hid_t plist_id  = H5Pcreate (H5P_DATASET_CREATE);
	H5Pset_chunk (plist_id, 1, dims);
	H5Pset_deflate (plist_id, 9);
	hid_t dataset_id = H5Dcreate2(file_id, dset_name.c_str(),datatype, dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	H5Dwrite (dataset_id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, names);	
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);	
	H5Pclose(plist_id);
	int row = 0;
	while (0 != (file_data = file->get (&errmsg))){
		dset_name = "/" + to_string(row);
		dataspace_id = H5Screate_simple(1, dims, NULL);
		if(dataspace_id < 0){
			cout << "Error in creating dataspace\n";
			return TCL_ERROR;
		}
		plist_id  = H5Pcreate (H5P_DATASET_CREATE);
		if(plist_id < 0){
			cout << "Error in creating plist\n";
			return TCL_ERROR;
		}		
		H5Pset_chunk (plist_id, 1, dims);
		H5Pset_deflate (plist_id, 9);		
		dataset_id = H5Dcreate2(file_id, dset_name.c_str(),datatype, dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
		if(dataset_id < 0){
			cout << "Error in creating dataset\n";
			return TCL_ERROR;
		}	
		H5Dwrite (dataset_id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, file_data);	
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		H5Pclose(plist_id);
		row++;
	}
	 H5Tclose(datatype);
	// H5Sclose(dataspace_id);
	 H5Fclose(file_id);
	 
	 
	 return TCL_OK;
}

extern "C" int convertHDF5Cmd(ClientData clientData, Tcl_Interp *interp,
                                         int argc,const char *argv[]){
	vector<string> args;
	
	for(int arg = 1; arg < argc; arg++){
		args.push_back(string(argv[arg]));
	}
	string str_arg;
	string input_file;
	string output_file;
	for(int arg = 0; arg < args.size(); arg++){
		str_arg = args[arg];
		transform(str_arg.begin(), str_arg.end(), str_arg.begin(), ::toupper);
		if((str_arg == "-I" || str_arg == "--I" || str_arg == "-INPUT" || str_arg == "--INPUT") && arg + 1 != args.size()){
			 input_file = args[++arg];
		 }else if((str_arg == "-O" || str_arg == "--O" || str_arg == "-OUTPUT" || str_arg == "--OUTPUT") && arg + 1 != args.size()) {
			 output_file =  args[++arg];
		 }else{
			cout << args[arg] << endl;
			return TCL_ERROR;
		}
	}
	cout << "here 1\n";
	return convert_file(interp, input_file.c_str(), output_file.c_str());										 
}
											 
											 

