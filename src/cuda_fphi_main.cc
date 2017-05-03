#include <cmath>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <random>
#include <cstdlib>
#include <iterator>
#include "solar.h"
#include <hdf5.h>
#include <chrono>
#include <Eigen/Dense>
using namespace std;

int generate_gpu_fphi_matrices(Tcl_Interp * interp, float ** trait_matrix, float ** cov_matrix, float ** aux, float ** eigenvectors,\
				vector<string> trait_names, int & n_subjects, int & n_covs);
extern int load_trait_matrix(Tcl_Interp * interp, Eigen::MatrixXd & Y, std::string  headers[], int dim);
static inline bool return_isspace (char  x){
	if(std::isspace(x))
		return true;
    else
	   return false;

}
static unsigned int * create_permutation_indices(int n_subjects, int n_permutations){
	unsigned int * permutation_indices = new unsigned int[n_subjects*n_permutations];
	for(int permutation = 0; permutation < n_permutations ; permutation++){
		iota (permutation_indices + permutation*n_subjects, permutation_indices + permutation*n_subjects + n_subjects, 0);
		random_shuffle (permutation_indices+ permutation*n_subjects, permutation_indices + permutation*n_subjects + n_subjects);	
	}		
	return permutation_indices;	
}

static unsigned int * create_permutation_matrix(unsigned int n_subjects, unsigned int n_permutations) {




  unsigned int * pmatrix = new unsigned int[n_subjects * n_permutations];

	for(unsigned int j = 0 ; j < n_permutations; j++){


		std::vector<unsigned int> randNums(n_subjects);

		std::iota (randNums.begin(), randNums.end(), 0);
		std::random_device rd;
		std::mt19937 mt(rd());

	    for (unsigned int n = n_subjects; n > 0; n--) {

	        std::uniform_int_distribution<unsigned int> dist(0.0,  n - 1);
	    	unsigned int Idx = dist(mt);
			pmatrix[(n - 1) + j*n_subjects] = randNums[Idx];
			randNums.erase(randNums.begin() + Idx);
	    }






	}



	return pmatrix;

}

static std::string convert_float_to_string(float value){
  std::stringstream ss;
  ss << value;
  std::string str_value(ss.str());

  return str_value;
}

static void write_connectivity_matrix(string file_out_name, string header_filename,float * h2 [], float * indicator [] ,float *  pvals[],bool  get_pval,size_t n_voxels){

	ifstream header_in(header_filename);

	vector<string> header_list;
	string current_label;
	while(header_in >> current_label){
		header_list.push_back(current_label);
	}
	header_in.close();
	string pval_name = ".pvalue_connectivity_matrix_" + file_out_name;
	string h2r_name = ".h2r_connectivity_matrix_" + file_out_name;

	string indicator_name = ".indicator_connectivity_matrix_" + file_out_name;

	ofstream h2r_conn_out (h2r_name.c_str());
	h2r_conn_out << ",";
	ofstream pval_conn_out;

	ofstream indicator_out (indicator_name.c_str());
	indicator_out << ",";

	if(get_pval){
		pval_conn_out.open(pval_name.c_str(), ofstream::out);
		pval_conn_out << ",";

	}

	for(vector<string>::iterator it = header_list.begin(); it != header_list.end(); it++){
		if((it + 1) == header_list.end()){
			h2r_conn_out << *it << "\n";
			indicator_out << *it << "\n";

			if(get_pval)
				pval_conn_out << *it << "\n";

		}else{
			h2r_conn_out << *it << ",";

			indicator_out << *it << ",";

			if(get_pval)
				pval_conn_out << *it << ",";
		}
	}


	for(size_t row = 0 ; row < n_voxels; row++){
		h2r_conn_out << header_list[row] << ",";
		indicator_out << header_list[row] << ",";
		if(get_pval)
			pval_conn_out << header_list[row] << ",";
		for(size_t col = 0 ; col < n_voxels ; col++){
			if((col + 1) == n_voxels){
				h2r_conn_out << h2[col][row]<< "\n";

				indicator_out << indicator[col][row] << "\n";

				if(get_pval)
					pval_conn_out << pvals[col][row] << "\n";
			}else{
				h2r_conn_out << h2[col][row] << ",";


				indicator_out << indicator[col][row] << ",";

				if(get_pval)
					pval_conn_out << pvals[col][row] << ",";
			}
		}
	}
	indicator_out.close();
	h2r_conn_out.close();

	if(get_pval)
		pval_conn_out.close();

}


static void write_column_line(string output_filename, string header_filename, size_t n_voxels){

	vector<string> header_list;
	string header;
	ifstream header_in(header_filename.c_str());

	while(header_in >> header)
		header_list.push_back(header);

	header_in.close();
	ofstream file_out(output_filename.c_str());
	file_out << ",";
	for(vector<string>::iterator it = header_list.begin(); it != header_list.end() ; it++){
		if(it+ 1 == header_list.end()){
			file_out << *it << "\n";
		}else{
			file_out << *it << ",";
		}

	}

	file_out.close();

}

static vector<string>  get_headers(string header_filename){
	vector<string> header_list;
	string header;
	ifstream header_in(header_filename.c_str());

	while(header_in >> header)
		header_list.push_back(header);

	header_in.close();

	return header_list;

}
static void write_row_data(string output_file, string headers[] ,  float * data [], size_t n_rows, size_t n_voxels){



	ofstream file_out(output_file.c_str(), fstream::out | fstream::app);
	for(size_t row = 0 ;row < n_rows ; row++){
		stringstream rows_out;
		rows_out << headers[row] << ",";
		for(size_t col = 0 ; col < n_voxels; col++){
			if(col + 1 == n_voxels)
				rows_out << data[row][col] << "\n";
			else
				rows_out << data[row][col] << ",";
		}

		file_out << rows_out.str();



	}

	file_out.close();

}

static void write_to_file(std::string file_name, std::string header_file, float * h2, float * indicator, float * SE, float * pval, bool get_pvalue, size_t n_voxels){

  std::ifstream file_in(header_file.c_str());
  std::vector< std::string > lines;

  if(file_in.is_open()){

	  for(size_t voxel = 0; voxel < n_voxels; voxel++){
		  std::string line;
		  if(file_in.eof()){
			  std::cout << "Warning the header file does not have enough trait names to correspond to the trait matrix column number.\n";


			  std::cout << "Proceeding anyway with trait number " << voxel + 1 << " out of " << n_voxels << " traits.\n";

			  line = ',' +  to_string(indicator[voxel]) + ',' + to_string(h2[voxel]) + ',' + to_string(SE[voxel]);
			  if(get_pvalue){
				  line +=  ',' + to_string(pval[voxel]) + '\n';
			  }else{
				  line += '\n';
			  }


		  }else{
			  file_in >> line;


		  	  line += ',' +  to_string(indicator[voxel]) + ',' + to_string(h2[voxel]) + ',' + to_string(SE[voxel]);
		  	  if(get_pvalue){
		  		  line +=  ',' + to_string(pval[voxel]) + '\n';
		  	  }else{
		  		  line += '\n';
		  	  }

		  }

		  lines.push_back(line);
	  }

  }else{
	  std::cout << "Warning header file was not found.  Proceeding anyway without trait names.\n";

	  for(size_t voxel = 0; voxel < n_voxels; voxel++){
		  std::string line;
		  line = "N/A,";
		  line +=  to_string(indicator[voxel]) + ',' + to_string(h2[voxel]);
		  if(get_pvalue){
			  line +=  ',' + to_string(pval[voxel]) + '\n';
		  }else{
			  line += '\n';
		  }
		  lines.push_back(line);
	  }
  }

  file_in.close();
  file_name = file_name + ".results.csv";
  std::ofstream file_out(file_name.c_str());
  file_out << "Trait,Indicator,H2r,SE";
  if(get_pvalue)
	  file_out << ",pvalue\n";
  else
	  file_out << "\n";

  for(std::vector< std::string >::iterator it = lines.begin() ; it != lines.end() ; it++){

	  file_out << *it;

  }

  file_out.close();

}
static void convert_matrix_to_pointer(float * ptr, Eigen::MatrixXd & matrix){
	
	for(int col = 0; col < matrix.cols(); col++) for(int row = 0; row < matrix.rows() ;row++) ptr[col*matrix.rows() + row] = matrix(row, col);
	
}
	

	
	
	
	
typedef struct{
	hid_t hdf5_file_out;

	hid_t matrix_dataset;

	hid_t matrix_filespace;
}hdf5_struct;
extern "C" std::vector<int> select_devices();
extern "C" std::vector<int> select_all_devices();
extern "C" int call_cudafphi(float * h2, float * indicator, float * SE, float * pvals, float * h_y,  const float * h_Z,  float * h_cov,
             const float * h_evectors, unsigned int * h_pmatrix, bool covariates, bool get_pval,
             size_t n_voxels, size_t n_subjects, size_t n_permutations, size_t n_covariates, std::vector<int> selected_devices);

extern "C" int call_cudafphi_connectivity(float * conn_h2, float * conn_indicator, float * conn_SE,
		 float * conn_pval, float * h_y,  const float * h_Z,  float * h_cov,
            const float * h_evectors,const unsigned int * h_pmatrix, bool covariates, bool get_pval, hdf5_struct * hdf5_data_h2, hdf5_struct * hdf5_data_indicator,
            hdf5_struct * hdf5_data_SE, size_t n_voxels, size_t n_subjects, size_t n_permutations, size_t n_covariates, std::vector<int> selected_devices);




void close_hdf5_ids(hdf5_struct * hdf5_data);

void initialize_hdf5_file_cuda_fphi(std::string headers [], const char * filename, hdf5_struct * hdf5_data, size_t row_begin, size_t row_end,
										size_t col_begin, size_t col_end);

void print_help(Tcl_Interp * interp){
	Solar_Eval(interp, "help gpu_fphi");	
}

extern "C" int gpufphiCmd (ClientData clientData, Tcl_Interp *interp,
		int argc, const char *argv[]){

	bool select_device = true;



	for (int index = 1 ; index < argc ; index++){
		string current_arg(argv[index]);
		string upper_arg = current_arg;
		transform(upper_arg.begin(), upper_arg.end(),upper_arg.begin(), ::toupper);
		if((upper_arg == "--HELP") || (upper_arg == "--H") || (upper_arg == "-HELP") || (upper_arg == "-H")){
			print_help(interp);
			return TCL_OK;
		}
	}
	std::vector<int> selected_devices;
	vector<string> arg_list;
	bool get_connectivity = false;

	for(int index = 1; index < argc; index++){
		arg_list.push_back(string(argv[index]));
	}
	for (vector<string>::iterator it = arg_list.begin(); it < arg_list.end(); it++){
		string current_arg(*it);
		string upper_arg = current_arg;
		transform(upper_arg.begin(), upper_arg.end(),upper_arg.begin(), ::toupper);
		if((upper_arg == "--ALL") || (upper_arg == "--A") || (upper_arg == "-ALL") || (upper_arg == "-A")){
			select_device = false;
			it = arg_list.erase(it);
			it--;
		}else if((upper_arg == "--CONN") || (upper_arg == "-CONN") || (upper_arg == "--C") || (upper_arg == "-C")){
			get_connectivity = true;
			it = arg_list.erase(it);
			it--;
		}
			
	}



	unsigned int n_permutations = 0;
	string file_out;
	string header_filename;
	for (int it = 0 ; it < arg_list.size(); it+=2){
		string current_arg = arg_list[it];
		string upper_arg = current_arg;
		transform(upper_arg.begin(), upper_arg.end(),upper_arg.begin(), ::toupper);
		if((upper_arg == "--NP"  || upper_arg == "-NP") && (it + 1 < arg_list.size())){
			n_permutations = strtol(arg_list[it + 1].c_str(), NULL, 10);
 			if(n_permutations <= 0){
 				cout << "Error in selecting the number of permutations.\n Must be greater than zero.\n";
 				return EXIT_FAILURE;
 			}
		}else if ((upper_arg == "--HEADER" || upper_arg == "-HEAD" || upper_arg == "--HEAD" || upper_arg == "-HEADER") &&
		    it + 1 < arg_list.size()){
			
			header_filename = arg_list[it + 1];	
			
		}else if (upper_arg == "--OUT" || upper_arg == "-OUT"  \
		|| upper_arg == "--O" || upper_arg == "-O" || upper_arg == "-OUTPUT" || upper_arg == "--OUTPUT"){
			
			file_out = arg_list[it + 1];	
			
		}else{

			cout << "Error in argument number " << it << ".\n";
			cout << "Argument: " << arg_list[it] << "\n";
			print_help(interp);
			return EXIT_FAILURE;

		}


	}
	
	if(select_device)
		selected_devices = select_devices();
	else
		selected_devices = select_all_devices();
	vector<string> header_list = get_headers(header_filename);
	if(header_list.size() == 0) {
		cout << "No traits were listed in the header.\n";
		return TCL_ERROR;
	}
	string headers[header_list.size()];
	for(vector<string>::iterator it = header_list.begin();
			it != header_list.end(); it++){
		size_t idx = distance(header_list.begin(), it );
		headers[idx] = *it;
	}
	
	if(header_filename == "" || file_out == ""){
		cout << "Missing a mandatory argument.\n";
		print_help(interp);
		return TCL_ERROR;
	}
	if(selected_devices.size() == 0){
		printf("No usable devices were selected.\n");
		return TCL_ERROR;
	}
	int  n_subjects;
	int  n_voxels = header_list.size();
	int n_covariates;
	Eigen::MatrixXd Y;
	//load_trait_matrix(interp, Y, headers, header_list.size());
	//Solar_Eval(interp, string("trait " + headers[0]).c_str()); 
	
	float * h_y, * cov, * aux, * evectors;
  std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
	
	int success = generate_gpu_fphi_matrices(interp, &h_y, &cov, &aux, &evectors,\
				header_list,  n_subjects,  n_covariates);
    if(success == TCL_ERROR) return success;
    bool covariates = false;
    if(n_covariates != 0){
		covariates = true;
		n_covariates++;
	}
	

	
 	bool get_pval = false;
 	unsigned int * pmatrix;
 	if(n_permutations > 0){
 		try{
 			pmatrix = create_permutation_indices(n_subjects, n_permutations);
 			get_pval = true;
		}catch(std::bad_alloc&){
			std::cout << "Could not allocate memory for permutation matrix.\n";
			delete [] h_y;
			delete [] evectors;
			delete [] aux;
			if(covariates)
				delete [] cov;
			return TCL_ERROR;
		}catch(...){
			std::cout << "Unknown failure in trying to create permutation matrix.\n";
			delete [] h_y;
			delete [] evectors;
			delete [] aux;			
			if(covariates)
				delete [] cov;
			return TCL_ERROR;
		}

 	}
	Y.resize(0, 0);
 	float * indicator;
 	float * pvals;
 	float * h2;
	float * SE;
 	try{
 	 	 indicator = new float[n_voxels];
		 if(get_pval)
			pvals = new float[n_voxels];
 	 	 h2 = new float[n_voxels];
 	 	 SE = new float[n_voxels];
 	}catch(std::bad_alloc&){
 		std::cout << "Error could not allocate memory for output arrays.\n";
	 	if(get_pval)
	 		delete [] pmatrix;
	 	delete [] h_y;
		delete [] evectors;
		delete [] aux;	
			
		if(covariates)
			delete [] cov;		 			
	 			
		return TCL_ERROR;
 	}
 
 	
 	if(!get_connectivity){
 		call_cudafphi(h2, indicator, SE, pvals, h_y,  aux,   cov,
 		              evectors, pmatrix,  covariates,  get_pval,
 		              n_voxels,  n_subjects,  n_permutations,  n_covariates,  selected_devices);

 		write_to_file(file_out,header_filename, h2, indicator,SE, pvals, get_pval, n_voxels);
		      
	}else{
 		hdf5_struct * hdf5_data_h2 = new hdf5_struct;
 		hdf5_struct * hdf5_data_indicator = new hdf5_struct;
 		hdf5_struct * hdf5_data_SE = new hdf5_struct;
 		hdf5_struct * hdf5_data_pvalues;
 		initialize_hdf5_file_cuda_fphi(header_list.data(), string(file_out + ".h2r.connectivity_matrix.h5").c_str(), hdf5_data_h2, 0, n_voxels,
 													0, n_voxels);
 		initialize_hdf5_file_cuda_fphi(header_list.data(), string(file_out + ".indicator.connectivity_matrix.h5").c_str(), hdf5_data_indicator, 0, n_voxels,
 													0, n_voxels);
 													
 		initialize_hdf5_file_cuda_fphi(header_list.data(), string(file_out + ".SE.connectivity_matrix.h5").c_str(), hdf5_data_SE, 0, n_voxels,
 													0, n_voxels);
		if(get_pval){	
			hdf5_data_pvalues = new hdf5_struct;												
			initialize_hdf5_file_cuda_fphi(header_list.data(), string(file_out + ".pvalues.connectivity_matrix.h5").c_str(), hdf5_data_pvalues, 0, n_voxels,
 													0, n_voxels);																						

		}

 		call_cudafphi_connectivity(h2, indicator, SE, pvals,  h_y, aux, cov, evectors, pmatrix, covariates,
 				get_pval,  hdf5_data_h2, hdf5_data_indicator, hdf5_data_SE, n_voxels, n_subjects, n_permutations, n_covariates, selected_devices);
 		close_hdf5_ids(hdf5_data_h2);
 		close_hdf5_ids(hdf5_data_indicator);
 		close_hdf5_ids(hdf5_data_SE);
 		if(get_pval){
			close_hdf5_ids(hdf5_data_pvalues);
			delete hdf5_data_pvalues;
		}
		delete hdf5_data_h2;
		delete hdf5_data_indicator;
		delete hdf5_data_SE;
	}
	
	delete [] h_y;
	delete [] evectors;
	delete [] aux;
	if(covariates)
		delete [] cov; 	
 	if(get_pval)
 		delete [] pmatrix;

 	delete [] h2;
 	delete [] SE;
 	if(get_pval)
 	 	delete [] pvals;
 	delete [] indicator;
 			
 	return TCL_OK;
}




