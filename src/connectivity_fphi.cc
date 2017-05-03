#include "solar.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <pthread.h>
#include <Eigen/Dense>
#include "hdf5.h"
#include <unordered_map>
#include <chrono>
#include <assert.h>
using namespace Eigen;
//using namespace std;
int generate_fphi_matrices(Tcl_Interp * interp, Eigen::MatrixXd  & l_trait_vector, Eigen::MatrixXd  & l_cov_matrix, Eigen::MatrixXd & l_aux, Eigen::MatrixXd & l_eigenvectors,\
				std::vector<std::string> trait_names, int n_traits);

extern MatrixXd* MathMatrixGet (char* mmid, const char* fname, Tcl_Interp* interp);
extern void calculate_h2r(MatrixXd Y, MatrixXd Z, MatrixXd ZTZI, float & h2r , float & score, float & SE);

static inline bool return_isspace (char  x){
	if(std::isspace(x))
		return true;
    else
	   return false;

}
static std::vector<std::string> parse_line(std::string currentLine, bool to_upper_case){

	currentLine.erase(std::remove_if(currentLine.begin(), currentLine.end(), return_isspace), currentLine.end());
	if(to_upper_case) transform(currentLine.begin(), currentLine.end(), currentLine.begin(), ::toupper);
	std::vector<std::string> output_vector;
	std::string current_str;
	char current_letter;
	for(auto it :currentLine){
		if(it == ','){
			output_vector.push_back(current_str);
			current_str = "";
		}else{
			current_str += it;
		}
	}

	output_vector.push_back(current_str);


	return output_vector;

}
static std::string parse_line_at_index(std::string currentLine, int index){

	currentLine.erase(std::remove_if(currentLine.begin(), currentLine.end(), return_isspace), currentLine.end());
	
	std::string current_str;
	char current_letter;
	int comma_count = 0;
	for(auto it :currentLine){
		
		
		if(it == ','){
			comma_count++;
			continue;
		}else if (comma_count == index){
			current_str += it;
		}else if (comma_count > index){
			return current_str;
		}
	}
	
	
	return current_str;


}




static void whithen_traits(MatrixXd & Y,  MatrixXd & Eigenvectors){

	MatrixXd Eigenvectors_T = Eigenvectors.transpose();
	for(size_t col = 0; col < Y.cols(); col++){
		VectorXd Y_col = Y.col(col);
		double mean = Y_col.mean();

		Y_col = (Y_col.array() - mean).matrix();
	
//		Y_col = (Y_col.array()/(norm/sqrt(Y_col.rows() - 1))).matrix();
		Y.col(col) = Eigenvectors_T*Y_col;

	}
	
}

static MatrixXd compute_hat(VectorXd  Y, MatrixXd X){
	MatrixXd cov(Y.rows(), X.cols() + 1);
	
	for(int col = 0 ; col < X.cols(); col++){
		cov.col(col) = X.col(col);
		
	}
	
	cov.col(cov.cols() - 1) = Y;
	
	return MatrixXd::Identity(Y.rows(), Y.rows()) - cov*(cov.transpose()*cov).inverse()*cov.transpose();
	
}
typedef struct{
 VectorXd Y;
 MatrixXd Z; 
	MatrixXd ZTZI; 
	MatrixXd hat; 
	float h2r; 
	float  score;
	float SE;
}fphi_data;	
static void * compute_residual_and_h2r(void * fphi_struct){
	
	fphi_data * data = (fphi_data*)fphi_struct;
	data->Y = data->hat*data->Y;
	if(data->Y.norm()/data->Y.rows() < 0.00001){
		data->h2r = 0.f;
		data->score = 0.f;
		data->SE = 0.f;
		return NULL;
	}

	calculate_h2r(data->Y, data->Z, data->ZTZI, data->h2r, data->score, data->SE);
}



extern "C" void write_columns_to_hdf5(float * data, hid_t & file_in, size_t n_cols, size_t row){
	
	
	hid_t dataset = H5Dopen(file_in, "/matrix", H5P_DEFAULT);
	
	hid_t filespace = H5Dget_space(dataset);
	hsize_t dimsext[2] = {1, n_cols};
	hsize_t offset[2] = {row, 0};	
	
    herr_t status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
				 dimsext, NULL);
				 
	assert(status >= 0);
	
	hid_t dataspace = H5Screate_simple(2, dimsext, NULL);	
	status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, dataspace, filespace,
		      H5P_DEFAULT, data);	
		      
		      
	assert(status >= 0);
	
	H5Dclose (dataset);
	
	H5Sclose (dataspace);
	H5Sclose (filespace);
	
	
}



static void compute_and_write_row_hdf5(float  * h2r, float * score, float * SE , MatrixXd & Y, MatrixXd & X, 
					 MatrixXd  & Z,  MatrixXd &  ZTZI,  size_t col_begin,
					size_t col_end, size_t row_iter, size_t row, hid_t & h2r_out , hid_t & score_out, hid_t & SE_out){
					
	size_t total_cols = col_end - col_begin;
	size_t col_iter = col_begin;
						 

	MatrixXd hat = compute_hat(Y.col(row_iter), X);
	
	fphi_data data[total_cols];
	for(size_t col = 0 ; col < total_cols; col++, col_iter++) {
		data[col].Y = Y.col(col_iter);
		data[col].Z = Z; 
		data[col].ZTZI = ZTZI; 				
		data[col].hat =  hat;  
		compute_residual_and_h2r(&data[col]);
		h2r[col] = data[col].h2r;
		score[col] = data[col].score;
		SE[col] = data[col].SE;			
	}
	
	write_columns_to_hdf5(h2r, h2r_out, total_cols, row);
	write_columns_to_hdf5(score, score_out, total_cols, row);
	write_columns_to_hdf5(SE, SE_out, total_cols, row);
		
				
}
static void write_csv_row(std::ofstream & file_out, float * data, size_t dim){
	
	for(size_t idx = 0; idx < dim; idx++) {		
		file_out << data[idx];		
		if(idx + 1 != dim) 
			file_out << ",";
		else
			file_out << "\n";
	}			
}
	
	
	
	
static void compute_and_write_row_csv(float * h2r, float  * score, float * SE, MatrixXd & Y, MatrixXd & X, 
					 MatrixXd  & Z,  MatrixXd &  ZTZI,  size_t col_begin,size_t col_end, size_t row, 
					 std::ofstream & h2r_out, std::ofstream & score_out, std::ofstream & SE_out){
						
	size_t total_cols = col_end - col_begin;
	size_t col_iter = col_begin;

	MatrixXd hat = compute_hat(Y.col(row), X);
	pthread_t threads[total_cols];
	fphi_data data[total_cols];
	for(size_t col = 0 ; col < total_cols; col++, col_iter++) {
		data[col].Y = Y.col(col_iter);
		data[col].Z = Z; 
		data[col].ZTZI = ZTZI; 				
		data[col].hat =  hat;  
			//	compute_residual_and_h2r(&data[col]);
		compute_residual_and_h2r(&data[col]);
		h2r[col] = data[col].h2r;
		score[col] = data[col].score;
		SE[col] = data[col].SE;		
	//	pthread_create(&threads[col], NULL, compute_residual_and_h2r, &data[col]);			
	}

	write_csv_row(h2r_out, h2r, total_cols);
	write_csv_row(score_out, score, total_cols);
	write_csv_row(SE_out, SE, total_cols);
		
				
}
typedef struct{
	MatrixXd  Y;
	MatrixXd  X;
	MatrixXd Z;
	MatrixXd ZTZI;
	size_t col_begin;
	size_t col_end;
	size_t row_iter;
	size_t row;
	
	float * h2r;
	float * score;
	float * SE;
	
	hid_t h2r_out;
	hid_t score_out;
	hid_t SE_out;
	
	
}fphi_row_data;
static void * compute_write_row_hdf5_pthread(void * pthread_data){
	
	fphi_row_data *row_data = (fphi_row_data*)pthread_data;
	
	size_t total_cols = row_data->col_end - row_data->col_begin;
	size_t col_iter = row_data->col_begin;
						 

	MatrixXd hat = compute_hat(row_data->Y.col(row_data->row_iter), row_data->X);
	fphi_data data;
	data.Z = row_data->Z; 
	data.ZTZI = row_data->ZTZI; 				
	data.hat =  hat;	
	for(size_t col = 0 ; col < total_cols; col++, col_iter++) {
		data.Y = row_data->Y.col(col_iter);
  
		compute_residual_and_h2r(&data);
		row_data->h2r[col] = data.h2r;
		row_data->score[col] = data.score;
		row_data->SE[col] = data.SE;			
	}
	
	write_columns_to_hdf5(row_data->h2r, row_data->h2r_out, total_cols, row_data->row);
	write_columns_to_hdf5(row_data->score,row_data->score_out, total_cols, row_data->row);
	write_columns_to_hdf5(row_data->SE,row_data->SE_out, total_cols, row_data->row);
			
	
	
}

static void calculate_connectivity_matrix(MatrixXd & Y, MatrixXd &  X, MatrixXd  & Z, std::ofstream & h2r_out_csv,
	std::ofstream & score_out_csv, std::ofstream & SE_out_csv,MatrixXd & eigenvectors, std::string headers [],size_t col_begin, size_t col_end, 
					size_t row_begin, size_t row_end, bool use_hdf5, hid_t & h2r_out, hid_t  & score_out, hid_t & SE_out){
						
	             
	whithen_traits(Y, eigenvectors);	

	size_t total_rows = row_end - row_begin;
	size_t row_iter = row_begin;
		
	MatrixXd  ZT = Z.transpose();
				
	MatrixXd  ZTZ = (ZT * Z);

	MatrixXd ZTZI = ZTZ.inverse();
	if(use_hdf5){
		int batch_size = 10;
		int iterations = ceil(float(total_rows)/batch_size);
		int free_traits = total_rows % batch_size;
		float ** h2r = new float*[batch_size];
		float ** score = new float*[batch_size];
		float ** SE = new float*[batch_size];
		pthread_t threads[batch_size];
		fphi_row_data data[batch_size];
		int row = 0;
		for(int batch = 0; batch < batch_size ; batch++){
			h2r[batch] = new float[1 + col_end - col_begin ];
		    score[batch] = new float[1 + col_end - col_begin ];
			SE[batch] = new float[1 + col_end - col_begin ];
		}			
		for(int iter = 0 ; iter < iterations ; iter++){
			
			for(int batch = 0; batch < batch_size && row < total_rows; batch++, row++, row_iter++){
				data[batch].Y = Y;
				data[batch].X = X;
				data[batch].ZTZI = ZTZI;
				data[batch].Z = Z;
				data[batch].col_begin = col_begin;
				data[batch].col_end = col_end;
				data[batch].row = row;
				data[batch].row_iter = row_iter;
				data[batch].h2r = h2r[batch];
				data[batch].score = score[batch];
				data[batch].SE = SE[batch];
				data[batch].h2r_out = h2r_out;
				data[batch].score_out = score_out;
				data[batch].SE_out = SE_out;
				
				pthread_create(&threads[batch], NULL, compute_write_row_hdf5_pthread, &data[batch]);
			}
			
			if(iter != iterations - 1){
				for(int batch = 0; batch<batch_size; batch++){
					pthread_join(threads[batch], NULL);
				}
			}else{
				for(int batch = 0; batch<free_traits; batch++){
					pthread_join(threads[batch], NULL);
				}	
			}				
		}
	//	for(size_t row = 0 ; row < total_rows; row++, row_iter++)\
	//	compute_and_write_row_hdf5(h2r, score, SE, Y, X, Z, ZTZI, col_begin, col_end, row_iter, row, h2r_out, score_out, SE_out);	
		for(int batch = 0; batch < batch_size ; batch++){
			delete [] h2r[batch];
			delete [] score[batch];
			delete [] SE[batch];
		}			
		delete [] h2r;
		delete [] score;
		delete [] SE;	
	}else{
		float * h2r = new float[1 + col_end - col_begin ];
		float  *score = new float [1 + col_end - col_begin ];
		float * SE = new float[1 + col_end - col_begin];		
		for(size_t row = 0; row < total_rows; row++ , row_iter++){
			h2r_out_csv << headers[row_iter] << ",";	
			score_out_csv << headers[row_iter] << ",";
			SE_out_csv << headers[row_iter] << ",";
			compute_and_write_row_csv( h2r ,  score, SE, Y,  X,  Z, ZTZI,   col_begin, col_end,  row_iter,  h2r_out_csv, score_out_csv, SE_out_csv);				
		}
		delete [] h2r;
		delete [] score;
		delete [] SE;										
								
	}
	
						
}


int initialize_hdf5_file(std::string headers [], const char * filename, hid_t  & hdf5_file, size_t row_begin, size_t row_end,
										size_t col_begin, size_t col_end){
											
	hsize_t col_dim = col_end - col_begin ;
	hsize_t row_dim = row_end - row_begin ;

	hsize_t label_row_dims[1] = {row_dim};
	hsize_t label_col_dims[1] = {col_dim};	
	
	std::string  row_headers[row_dim];
	size_t row_iter = row_begin;
	for(size_t idx = 0; idx < row_dim ; idx++) row_headers[idx] = headers[row_iter++].c_str();
	std::string col_headers[col_dim];
	size_t col_iter = col_begin;
	for(size_t idx = 0; idx < col_dim ; idx++) col_headers[idx] = headers[col_iter++].c_str();	
	
	
    hid_t datatype = H5Tcopy (H5T_C_S1);
    assert (datatype >= 0);
    
    herr_t status = H5Tset_size(datatype, H5T_VARIABLE);
    
    assert (status >= 0);
    
    status = H5Tset_strpad(datatype, H5T_STR_NULLTERM);
    
    hid_t row_dataspace = H5Screate_simple (1, label_row_dims, NULL);

  //  status = H5Tset_size (datatype,H5T_VARIABLE);
  //  assert (status >= 0);
	
	hid_t row_plist_id  = H5Pcreate (H5P_DATASET_CREATE);

    status = H5Pset_chunk (row_plist_id, 1, label_row_dims);
    assert (status >= 0);

    status = H5Pset_deflate (row_plist_id, 9); 		
    assert (status >= 0);
    
    hid_t row_dataset_id = H5Dcreate2 (hdf5_file, "/row_headers", datatype, 
                            row_dataspace, H5P_DEFAULT, row_plist_id, H5P_DEFAULT);
                            
                            
    status = H5Dwrite (row_dataset_id, datatype, row_dataspace, H5S_ALL, H5P_DEFAULT, row_headers);
    assert(status >= 0);
     
    H5Sclose (row_dataspace);
    H5Dclose (row_dataset_id);
    H5Pclose (row_plist_id);
     
     
	hid_t col_dataspace = H5Screate_simple (1, label_col_dims, NULL);
		
	
	hid_t col_plist_id  = H5Pcreate (H5P_DATASET_CREATE);

    status = H5Pset_chunk (col_plist_id, 1, label_col_dims);
    assert (status >= 0);

    status = H5Pset_deflate (col_plist_id, 9); 		
    assert (status >= 0);
    
    hid_t col_dataset_id = H5Dcreate2 (hdf5_file, "/column_headers", datatype, 
                            col_dataspace, H5P_DEFAULT, col_plist_id, H5P_DEFAULT);
                            
                            
    status = H5Dwrite (col_dataset_id, datatype, col_dataspace, H5S_ALL, H5P_DEFAULT, col_headers);
    assert(status >= 0);
     
    H5Sclose (col_dataspace);
    H5Dclose (col_dataset_id);
    H5Pclose (col_plist_id);  
    H5Tclose (datatype);
     
	hsize_t matrix_dims[2] = {row_dim, col_dim};
	hsize_t matrix_chunk_dims[2] = {1, col_dim};  
     
    hid_t matrix_space =  H5Screate_simple (2, matrix_dims, NULL);
	hid_t matrix_plist_id  = H5Pcreate (H5P_DATASET_CREATE);

	status = H5Pset_chunk (matrix_plist_id, 2, matrix_chunk_dims);
    assert (status >= 0);

    status = H5Pset_deflate (matrix_plist_id, 9); 		
    assert (status >= 0); 
    const float fill_value = 0.f;
    status = H5Pset_fill_value(matrix_plist_id, H5T_NATIVE_FLOAT, &fill_value);
	assert (status >= 0); 
	 
    hid_t matrix_dataset_id = H5Dcreate2 (hdf5_file, "/matrix", H5T_NATIVE_FLOAT, 
                            matrix_space, H5P_DEFAULT, matrix_plist_id, H5P_DEFAULT);
                            
    H5Sclose(matrix_space);
    H5Dclose(matrix_dataset_id);
    H5Pclose(matrix_plist_id);
     
     
    return TCL_OK;
}

std::vector<std::string> read_headers(const char * filename){
	
	
	std::string current_header;
	std::ifstream file_in(filename);
	std::vector<std::string> buffer;
	while(file_in >> current_header) buffer.push_back(current_header);

	file_in.close();
	
	return buffer;
	
	
}
static void print_conn_fphi_help(Tcl_Interp * interp){
		Solar_Eval(interp, "help connectivty_fphi");
}

static void initialize_matrix_csv(std::ofstream &file_out, std::string headers[],
									size_t col_begin, size_t col_end){
	file_out << ",";
	for(size_t idx = col_begin ; idx < col_end; idx++){
		file_out << headers[idx];
		if(idx != col_end)
			file_out << ",";
		else
			file_out << "\n";
	}
			
}
int load_trait_matrix(Tcl_Interp * interp, MatrixXd & Y, std::string  headers[], int dim){
	std::string phenotype_file = std::string(Phenotypes::filenames());	
	
	std::ifstream pheno_in(phenotype_file.c_str());
	std::string line;
	std::getline(pheno_in, line);
	
	std::vector<std::string> title_line = parse_line(line, true);
	int id_idx = std::distance(title_line.begin(), std::find(title_line.begin(), title_line.end(), "ID"));
	int header_indices[dim];
	
	for(int idx =0; idx < dim;idx++){
		std::transform(headers[idx].begin(), headers[idx].end(), headers[idx].begin(), ::toupper);
		header_indices[idx] = std::distance(title_line.begin(), std::find(title_line.begin(), title_line.end(), headers[idx]));
	}
	std::unordered_map<std::string, std::vector<double> > id_map;
	std::vector<std::string> ids;
	std::vector<std::string> parsed_line;
	while(getline(pheno_in, line)){
		parsed_line = parse_line(line, false);
		
		std::vector<double> traits_at_row(dim);
		bool skip_row = false;
		for(int idx =0 ; idx < dim; idx++){
			std::string value = parsed_line[header_indices[idx]];
			if(value == ""){
				skip_row = true;
				break;
			}
			
			traits_at_row[idx] = atof(value.c_str());
			
		}
		
		if(skip_row) continue;
		
		id_map.insert(std::pair<std::string, std::vector<double> >(parsed_line[id_idx], traits_at_row));
	}
	
	pheno_in.close();
	
	
	std::ifstream pedindex_in("pedindex.csv");
	
	if(pedindex_in.is_open() == false){
		Solar_Eval(interp, "ped2csv pedindex.out pedindex.csv");
		pedindex_in.open("pedindex.csv", std::ifstream::in);
	}
	
	getline(pedindex_in, line);
	title_line = parse_line(line, true);
	id_idx = std::distance(title_line.begin(), std::find(title_line.begin(), title_line.end(), "ID"));
	std::vector< double * > all_rows; 
	std::string id; 
	std::unordered_map<std::string, std::vector<double> >::iterator map_iter;
	std::vector<double> current_row;
	while(std::getline(pedindex_in, line)){
		id = parse_line_at_index(line, id_idx);
		
		map_iter = id_map.find(id);
		
	
		
		if(map_iter != id_map.end()) all_rows.push_back(map_iter->second.data());
			
	}
	
	pedindex_in.close();
	MatrixXd temp_Y(all_rows.size(), dim);
	for(int row_idx = 0 ; row_idx < all_rows.size(); row_idx++) \
		for(int col_idx =0 ; col_idx < dim; col_idx++) temp_Y(row_idx, col_idx) = all_rows[row_idx][col_idx];
	

	Y = temp_Y;
	return TCL_OK;
}	
	
extern "C" int Runconnfphicmd(ClientData clientData, Tcl_Interp *interp,
                 int argc,const char *argv[]){
					 
					 
	Eigen::initParallel();				 
                 
					 
	std::vector<std::string> args;
	
	for(int arg = 1 ; arg < argc ; arg++) args.push_back(std::string(argv[arg]));
	MatrixXd  Z;
	MatrixXd  Y;
	MatrixXd  X;
	MatrixXd  eigenvectors;	
	
	const char * header_filename;
	const char * output_filename;
	bool use_hdf5 = true;
	size_t row_begin = 0, row_end = 0, col_begin = 0, col_end = 0;
	for(std::vector<std::string>::iterator it = args.begin() ;  it != args.end(); it++){
		transform((*it).begin(), (*it).end(), (*it).begin(), ::toupper);
		if((*it).compare("-H") == 0 || (*it).compare("--H") == 0 ||
				(*it).compare("-HELP") == 0 || (*it).compare("--HELP") == 0 ){
			print_conn_fphi_help(interp);
			return TCL_OK;
		}else if(((*it).compare("-ROWS") == 0 ||
				(*it).compare("--ROWS") == 0) && (it + 2 != args.end())){
			row_begin = std::stoi((*(++it)).c_str());
			row_end =  std::stoi((*(++it)).c_str());	
		}else if(((*it).compare("-COLS") == 0 ||
				(*it).compare("--COLS") == 0) && (it + 2 != args.end())){
			col_begin = std::stoi((*(++it)).c_str());
			col_end =  std::stoi((*(++it)).c_str());	
		}else if(((*it).compare("-O") == 0 ||(*it).compare("-OUT") == 0 ||
				(*it).compare("--O") == 0 || (*it).compare("--OUT") == 0 ) 
				&& (it + 1 != args.end())){
			output_filename = (*(++it)).c_str();		
		}else if(((*it).compare("-HEADER") == 0 ||(*it).compare("-HEAD") == 0 ||
				(*it).compare("--HEADER") == 0 || (*it).compare("--HEAD") == 0 ) 
				&& (it + 1 != args.end())){
			header_filename = (*(++it)).c_str();		
		}else if((*it).compare("-CSV") == 0 ||
				(*it).compare("--CSV") == 0  ){
			use_hdf5 = false;
		}else{
			fprintf(stderr, "An invalid argument was entered.\n");
			print_conn_fphi_help(interp);
			return TCL_ERROR;
		}
	}
	
	std::vector< std::string > buffer = read_headers(header_filename);
	std::string headers[buffer.size()];
	for(size_t idx = 0 ; idx < buffer.size(); idx++) headers[idx] = buffer[idx];	
	
	//load_trait_matrix(interp, Y, headers, buffer.size());
	//Solar_Eval(interp, std::string("trait " + headers[0]).c_str());
	//Eigen::VectorXd dummy;
	generate_fphi_matrices(interp, Y, X, Z, eigenvectors, buffer, buffer.size());	 

	
	if(row_end == 0) row_end = Y.cols();
	if(col_end == 0) col_end = Y.cols();	
	
	
	
	if(row_end < row_begin || col_end < col_begin || 
	col_begin < 0 || col_end < 0 || row_begin < 0 || row_end < 0){
		fprintf(stderr, "The region selected was invalid.\n");
		print_conn_fphi_help(interp);
		return TCL_ERROR;
	}
	
	if(Y.rows() != Z.rows()){
		fprintf(stderr, "The number of subjects of the trait matrix (%i) is not equal to that of the eigenvectors (%i).\n", Y.rows(), Z.rows());
		return TCL_ERROR;
	}

	
											 

	
	std::string h2r_filename;
	std::string indicator_filename; 	
	std::string SE_filename;			 							 

	hid_t hdf5_h2r_output;
	hid_t  hdf5_indicator_output;
	hid_t hdf5_SE_output;
	
	
	std::ofstream h2r_out_csv;
	std::ofstream score_out_csv;
	std::ofstream SE_out_csv;
	if(use_hdf5){		
		h2r_filename = std::string(output_filename) + ".h2r_connectivty_matrix.h5";
		 
		indicator_filename = std::string(output_filename) + ".indicator_connectivity_matrix.h5";
		
		SE_filename = std::string(output_filename) + ".SE.connectivity_matrix.h5";
		
		hdf5_h2r_output = H5Fcreate(h2r_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		
		
		hdf5_indicator_output = H5Fcreate(indicator_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		
		hdf5_SE_output = H5Fcreate(SE_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		
		 
		if(initialize_hdf5_file(headers, h2r_filename.c_str(), hdf5_h2r_output, row_begin, row_end, \
										 col_begin,  col_end) == TCL_ERROR) return TCL_ERROR;	
		if(initialize_hdf5_file(headers, indicator_filename.c_str(), hdf5_indicator_output, row_begin, row_end, \
									 col_begin,  col_end) == TCL_ERROR) return TCL_ERROR;	
									 
		if(initialize_hdf5_file(headers, SE_filename.c_str(), hdf5_SE_output, row_begin, row_end, \
					 col_begin,  col_end) == TCL_ERROR) return TCL_ERROR;							 			
	}else{
		h2r_filename = std::string(output_filename) + ".h2r.connectivty_matrix.csv";
		 
		indicator_filename = std::string(output_filename) + ".indicator.connectivity_matrix.csv";	
		SE_filename = std::string(output_filename) + ".SE.connectivity_matrix.csv";
		h2r_out_csv.open(h2r_filename.c_str(), std::fstream::out);
		
		score_out_csv.open(indicator_filename.c_str(), std::fstream::out);	
		initialize_matrix_csv(h2r_out_csv,  headers,  col_begin, col_end);
		initialize_matrix_csv(score_out_csv,  headers, col_begin, col_end);
		initialize_matrix_csv(SE_out_csv,  headers, col_begin, col_end);
	}									 
	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();				 				 
	calculate_connectivity_matrix( Y,  X,  Z,  h2r_out_csv, score_out_csv, SE_out_csv, eigenvectors,  headers, col_begin,  col_end, 
					 row_begin, row_end,  use_hdf5,  hdf5_h2r_output,   hdf5_indicator_output, hdf5_SE_output);
					 
	auto elapsed = std::chrono::high_resolution_clock::now() - start;
	std::chrono::duration<double> seconds = std::chrono::duration_cast<std::chrono::duration<double>>(elapsed);					 
	printf("%f \n", seconds.count());
	if(use_hdf5){
		H5Fclose( hdf5_h2r_output);
		H5Fclose( hdf5_indicator_output);
		H5Fclose( hdf5_SE_output);
	}else{
		h2r_out_csv.close();
		score_out_csv.close();
		SE_out_csv.close();
	}
		
					 

	
	return TCL_OK;
						 
}
	
	
	
	

