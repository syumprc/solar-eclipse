#include<cmath>
#include<stdio.h>
#include<vector>
#include<algorithm>
#include<sstream>
#include<fstream>
#include<iostream>
#include<string>
#include <time.h>
#include<random>
#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>

using namespace std;
static bool return_isspace (char  x){
	if(std::isspace(x)){
		return true;
	}
	return false;
}




static std::vector<float> parse_line(std::string currentLine){

	currentLine.erase(std::remove_if(currentLine.begin(), currentLine.end(), return_isspace), currentLine.end());



	std::vector<float> output_vector;
	std::string current_str;
	char current_letter;

	for(std::string::iterator it = currentLine.begin(); it != currentLine.end(); it++){
		current_letter = *it;
		if(current_letter == ','){
			output_vector.push_back(std::strtod(current_str.c_str(), NULL));
			current_str.clear();
		}else{
			current_str += current_letter;
		}
	}
	output_vector.push_back(std::strtod(current_str.c_str(), NULL));


	return output_vector;

}



static float * csvread_float(const char * filename, size_t &rows,
	    size_t &columns) {


	std::ifstream inputData(filename);
;
	if(inputData.is_open() == false){
		std::cerr << "Error in opening file: " << filename << "\n";
		return NULL;
	}
	std::string Line;

	std::vector< std::vector<float> > data;

	while (std::getline(inputData, Line)) {

		data.push_back(parse_line(Line));

	}

	rows = data.size();

	columns = data[0].size();

	inputData.close();

	float * data_out = new float [rows*columns];


	size_t colIdx = 0;
	size_t rowIdx = 0;
	for(std::vector< std::vector<float> >::iterator iter = data.begin(); iter != data.end() ; iter++, rowIdx++){
		colIdx = 0;
		for(std::vector<float>::iterator current_value = iter->begin() ; current_value != iter->end(); current_value++, colIdx++){
			data_out[colIdx*rows + rowIdx] = *current_value;
		}
	}


	return data_out;

}



static unsigned int * create_permutation_matrix(unsigned int n_subjects, unsigned int n_permutations) {




  unsigned int * pmatrix = new unsigned int[n_subjects * n_permutations];
  srand(0);
	for(unsigned int j = 0 ; j < n_permutations; j++){


		std::vector<unsigned int> randNums;

		for (int r = 0; r < n_subjects; r++) {
			randNums.push_back(r);
		}
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
	string pval_name = "pvalue_connectivity_matrix_" + file_out_name;
	string h2r_name = "h2r_connectivity_matrix_" + file_out_name;

	string indicator_name = "indicator_connectivity_matrix_" + file_out_name;

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

static void get_headers(string header_filename, string headers [], int node_number, size_t qsub_shared_batch_size, size_t qsub_batch_size){
	vector<string> header_list;
	string header;
	ifstream header_in(header_filename.c_str());

	while(header_in >> header)
		header_list.push_back(header);

	header_in.close();

	for(vector<string>::iterator it = header_list.begin() +  node_number * qsub_shared_batch_size ;
			it != header_list.begin() +  node_number * qsub_shared_batch_size + qsub_batch_size; it++){
		size_t idx = distance(header_list.begin() + node_number*qsub_shared_batch_size, it );
		headers[idx] = *it;
	}

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

static void write_to_file(std::string file_name, std::string header_file, float * h2, float * indicator, float * pval, bool get_pvalue, size_t n_voxels){

  std::ifstream file_in(header_file.c_str());
  std::vector< std::string > lines;

  if(file_in.is_open()){

	  for(size_t voxel = 0; voxel < n_voxels; voxel++){
		  std::string line;
		  if(file_in.eof()){
			  std::cout << "Warning the header file does not have enough trait names to correspond to the trait matrix column number.\n";


			  std::cout << "Proceeding anyway with trait number " << voxel + 1 << " out of " << n_voxels << " traits.\n";

			  line = ',' +  convert_float_to_string(indicator[voxel]) + ',' + convert_float_to_string(h2[voxel]);
			  if(get_pvalue){
				  line +=  ',' + convert_float_to_string(pval[voxel]) + '\n';
			  }else{
				  line += '\n';
			  }


		  }else{
			  file_in >> line;


		  	  line += ',' +  convert_float_to_string(indicator[voxel]) + ',' + convert_float_to_string(h2[voxel]);
		  	  if(get_pvalue){
		  		  line +=  ',' + convert_float_to_string(pval[voxel]) + '\n';
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
		  line +=  convert_float_to_string(indicator[voxel]) + ',' + convert_float_to_string(h2[voxel]);
		  if(get_pvalue){
			  line +=  ',' + convert_float_to_string(pval[voxel]) + '\n';
		  }else{
			  line += '\n';
		  }
		  lines.push_back(line);
	  }
  }

  file_in.close();
  file_name = file_name + "_results.csv";
  std::ofstream file_out(file_name.c_str());
  file_out << "Trait,Indicator,H2r";
  if(get_pvalue)
	  file_out << ",pvalue\n";
  else
	  file_out << "\n";

  for(std::vector< std::string >::iterator it = lines.begin() ; it != lines.end() ; it++){

	  file_out << *it;

  }

  file_out.close();

}

extern "C" std::vector<int> select_devices();
extern "C" int call_cudafphi(float * h2, float * indicator, float * pvals, float * conn_h2, float * conn_indicator,
		 float * conn_pval, float * h_y,  const float * h_Z, const float * h_cov,
            const float * h_evectors,const unsigned int * h_pmatrix, bool covariates, bool get_pval, bool get_connectivity,
            size_t n_voxels, size_t n_subjects, size_t n_permutations, size_t n_covariates, std::vector<int> selected_devices);

extern "C" int call_cudafphi_qsub(float * conn_h2, float * conn_indicator ,
		 float * conn_pval , float * h_y,  const float * h_Z,  float * h_cov,
            const float * h_evectors,const unsigned int * h_pmatrix, bool covariates, bool get_pval, bool get_connectivity, int node_number,
            size_t n_voxels, size_t n_subjects, size_t n_permutations, size_t qsub_batch_size, size_t qsub_shared_batch_size,
            size_t n_covariates, std::vector<int> selected_devices);
extern "C" std::vector<int> select_all_devices();


extern void write_file(const char * out_filename, string labels [], float **  data , size_t dim);
int write_to_hdf5_qsub(const char * out_filename, float * data, size_t qsub_batch_size, size_t qsub_shared_batch_size, size_t node_number, size_t dim);

void write_data_to_hdf5(const char * out_filename, string  labels [], float *  buffer , size_t dim);

int write_to_hdf5_qsub(string out_base_name,  float * h2r , float * indicator , float * pvals ,size_t qsub_batch_size, size_t qsub_shared_batch_size, size_t node_number, size_t n_voxels, bool use_pval){
	string h2r_filename = out_base_name + "_h2r_connectivity_matrix.h5";
	string indicator_filename = out_base_name + "_indicator_connectivity_matrix.h5";

	int error;

	error = write_to_hdf5_qsub(h2r_filename.c_str(), h2r,  qsub_batch_size,  qsub_shared_batch_size,  node_number,  n_voxels);

	if(error == 0)
		return error;

	error = write_to_hdf5_qsub(indicator_filename.c_str(), indicator,   qsub_batch_size,  qsub_shared_batch_size,  node_number, n_voxels);
	if(error == 0)
		return error;

	if(use_pval){
		string pval_filename = out_base_name + "_pvalues_connectivity_matrix.h5";
		error = write_to_hdf5_qsub(pval_filename.c_str(),pvals,   qsub_batch_size,  qsub_shared_batch_size,  node_number,  n_voxels);
		if(error == 0)
			return error;
	}


	return 1;

}










void print_help(){
	printf("cuda_fphi calculates the heritability and indicator values of a set of traits that share an auxiliary and eigenvector matrix.\n");
	printf("Optionally the pvalue can be calculated given a number of permutations and covariates can be included given a hat matrix.\n");
	printf("Where hat = I - X * (XT * X)^-1 * XT\n");
	printf("cuda_fphi -Y <trait matrix file name>  -aux <auxiliary matrix filename> -eigen <eigenvector matrix (non-transposed) filename>  -o <output file name>  -header <trait header file name> optional args: <-conn connectivity matrix switch> <-all switch> -np <n_permutations> -cov <covariate matrix file name>\n");
	printf("To combine cuda_fphi and qsub simply add -qsub <job id> <total jobs>\n");
	printf("An error log called 'cuda_fphi.log' will can create with the --debug option\n");
}

int main (int argc, const char *argv[]){

	bool select_device = true;
	bool use_qsub = false;
	close(STDERR_FILENO);
	int error_file;

	time_t rawtime;
	struct tm * timeinfo;
	bool debug = false;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	for (int index = 1 ; index < argc ; index++){
		string current_arg(argv[index]);
		string upper_arg = current_arg;
		transform(upper_arg.begin(), upper_arg.end(),upper_arg.begin(), ::toupper);
		if((upper_arg == "--HELP") || (upper_arg == "--H") || (upper_arg == "-HELP") || (upper_arg == "-H")){
			print_help();
			return EXIT_SUCCESS;
		}
	}

	std::vector<int> selected_devices;
	vector<string> arg_list;
	bool get_connectivity = false;

	size_t node_number;
	size_t total_nodes;

	for(int index = 1; index < argc; index++){
		arg_list.push_back(string(argv[index]));
	}
	for (vector<string>::iterator it = arg_list.begin(); it < arg_list.end(); it++){
		string current_arg(*it);
		string upper_arg = current_arg;
		transform(upper_arg.begin(), upper_arg.end(),upper_arg.begin(), ::toupper);
		if((upper_arg == "--ALL") || (upper_arg == "--A") || (upper_arg == "-ALL") || (upper_arg == "-A")){
			selected_devices = select_all_devices();
			select_device = false;
			it = arg_list.erase(it);
			it--;
		}else if((upper_arg == "--CONN") || (upper_arg == "-CONN") || (upper_arg == "--C") || (upper_arg == "-C")){
			get_connectivity = true;
			it = arg_list.erase(it);
			it--;
		}else if((upper_arg == "--QSUB" || upper_arg == "-QSUB") &&  (it + 2 != arg_list.end())){
			use_qsub = true;
			
			int arg_index =  distance(arg_list.begin(), it);
			node_number = strtol(arg_list[arg_index + 1].c_str(), NULL, 10);
			total_nodes = strtol(arg_list[arg_index + 2].c_str(), NULL, 10);
			it = arg_list.erase(it, it + 3);
			it--;
		}else if((upper_arg == "--DEBUG" || upper_arg == "-DEBUG" || upper_arg == "-D" || upper_arg == "--D") &&  (it + 2 != arg_list.end())){
			it = arg_list.erase(it);
			debug = true;
			error_file = open( "cuda_fphi.log", O_CREAT | O_WRONLY, 0644 );

			dup2( error_file, STDERR_FILENO);

			fprintf(stderr, "CUDA FPHI Error Log\n");
			fprintf (stderr, "Current local time and date: %s", asctime (timeinfo) );
			it--;
		}
	}


	if(select_device)
		selected_devices = select_devices();



	string Y_file;
	string evectors_file;
	string aux_file;
	string cov_file;
	unsigned int n_permutations = 0;
	string file_out;
	string header_filename;
	for (int it = 0 ; it < arg_list.size(); it+=2){
		string current_arg = arg_list[it];
		string upper_arg = current_arg;
		transform(upper_arg.begin(), upper_arg.end(),upper_arg.begin(), ::toupper);
		if((upper_arg == "--Y"  || upper_arg == "-Y") && (it + 1 < arg_list.size())){
			Y_file = arg_list[it + 1];
		}else if((upper_arg == "--AUX"  || upper_arg == "-AUX") && (it + 1 < arg_list.size())){
			aux_file  = arg_list[it + 1];
		}else if((upper_arg == "--EIGEN"  || upper_arg == "-EIGEN" || upper_arg == "-E" || upper_arg == "--E") && (it + 1 < arg_list.size())){
			evectors_file = arg_list[it + 1];
		}else if((upper_arg == "--O"  || upper_arg == "-O" || upper_arg == "--OUT" || upper_arg == "-OUT" || upper_arg == "--OUTPUT" || upper_arg == "-OUTPUT") && (it + 1 < arg_list.size())){
			file_out = arg_list[it + 1];
		}else if((upper_arg == "--HEAD"  || upper_arg == "-HEAD" || upper_arg == "--HEADER" || upper_arg == "-HEADER") && (it + 1 < arg_list.size())){
			header_filename = arg_list[it + 1];
		}else if((upper_arg == "--NP"  || upper_arg == "-NP") && (it + 1 < arg_list.size())){
			n_permutations = strtol(arg_list[it + 1].c_str(), NULL, 10);
 			if(n_permutations <= 0){
 				cout << "Error in selecting the number of permutations.\n Must be greater than zero.\n";
 				return EXIT_FAILURE;
 			}
		}else if((upper_arg == "--COV"  || upper_arg == "-COV" || upper_arg == "-COVARIATES" || upper_arg == "--COVARIATES") && (it + 1 < arg_list.size())){
			cov_file = arg_list[it + 1];
		}else{

			cout << "Error in argument number " << it << ".\n";
			cout << "Argument: " << arg_list[it] << "\n";
			print_help();
			return EXIT_FAILURE;

		}


	}


	if(aux_file == "" || Y_file == "" || evectors_file == "" || file_out == "" || header_filename == ""){

		cout << "Missing a mandatory argument.\n";
		print_help();
		return EXIT_FAILURE;
	}


	if(selected_devices.size() == 0){
		printf("No usable devices were selected.\n");
		return EXIT_FAILURE;
	}



	size_t  n_subjects;
	size_t  n_voxels;
	float * h_y;

	try{
		h_y = csvread_float(Y_file.c_str(), n_subjects, n_voxels);
	}catch(std::bad_alloc&){
		std::cout << "Failed to allocate the memory needed for the trait matrix.\n";
		return EXIT_FAILURE;
	}catch(...){
		std::cout << "Unknown failure in trying to load the trait matrix.\n";
		return EXIT_FAILURE;
	}


 	if(h_y == NULL){
 		return EXIT_FAILURE;
 	}

 	size_t n_subjects_1;
 	size_t n_subjects_2;
 	float * evectors;


	try{

		evectors = csvread_float(evectors_file.c_str(), n_subjects_1, n_subjects_2);

	}catch(std::bad_alloc&){
		std::cout << "Failed to allocate the memory needed for the eigenvector matrix.\n";
		delete [] h_y;
		return EXIT_FAILURE;
	}catch(...){
		std::cout << "Unknown failure in trying to load the eigenvector matrix.\n";
		delete [] h_y;
		return EXIT_FAILURE;
	}

 	if(evectors == NULL){
 		delete [] h_y;
 		return EXIT_FAILURE;
 	}

 	if(n_subjects_1 != n_subjects_2){
 		printf("Row size and column size are not equal for evector matrix.\n There are %i rows and %i columns.\n",
			  n_subjects_1, n_subjects_2);
 		delete [] h_y;
 		delete [] evectors;
 		return EXIT_FAILURE;
 	}

 	if(n_subjects_1 != n_subjects){
 		printf("Rows of trait matrix and eigenvector matrix are not equal.\n %i rows in the trait matrix and %i rows and columns in the eigenvector matrix.\n",
			  n_subjects, n_subjects_1);
 		delete [] h_y;
 		delete [] evectors;
 		return EXIT_FAILURE;
 	}

 	size_t  dummy;
 	float * cov;
 	float *  aux;


	try{
		aux = csvread_float(aux_file.c_str(), n_subjects_1, dummy);
	}catch(std::bad_alloc&){
		std::cout << "Failed to allocate the memory needed for the auxiliary matrix.\n";
		delete [] evectors;
		delete [] h_y;
		return EXIT_FAILURE;
	}catch(...){
		std::cout << "Unknown failure in trying to load the auxiliary matrix.\n";
		delete [] evectors;
		delete [] h_y;
		return EXIT_FAILURE;
	}


 	if(aux == NULL){
 		delete [] evectors;
 		delete [] h_y;
 		return EXIT_FAILURE;
 	}

 	if(dummy != 2){
 		printf("Auxiliary matrix requires only two columns.  %i columns were detected.\n", dummy);
 		delete [] aux;
 		delete [] evectors;
 		delete [] h_y;
 		return EXIT_FAILURE;
 	}

 	if(n_subjects_1 != n_subjects){
 		printf("Auxiliary matrix subject number is %i while the trait matrix subject number is %i.\n",n_subjects_1, n_subjects);
 		delete [] aux;
 		delete [] evectors;
 		delete [] h_y;
 		return EXIT_FAILURE;
 	}


 	bool covariates = false;
 	bool get_pval = false;
 	unsigned int * pmatrix;

 	size_t  n_covariates;

 	if(n_permutations > 0){
 		try{
 			pmatrix = create_permutation_matrix(n_subjects, n_permutations);
 			get_pval = true;
		}catch(std::bad_alloc&){
			std::cout << "Could not allocate memory for permutation matrix.\n";
	 		delete [] aux;
	 		delete [] evectors;
	 		delete [] h_y;
			return EXIT_FAILURE;
		}catch(...){
			std::cout << "Unknown failure in trying to create permutation matrix.\n";
	 		delete [] aux;
	 		delete [] evectors;
	 		delete [] h_y;
			return EXIT_FAILURE;
		}

 	}

 	if(cov_file != ""){
 		try{
 			cov = csvread_float(cov_file.c_str(), n_subjects_1, n_covariates);
 			covariates = true;
		}catch(std::bad_alloc&){
			std::cout << "Could not allocate memory for covariate matrix matrix.\n";
	 		delete [] aux;
	 		delete [] evectors;
	 		delete [] h_y;
	 		if(get_pval)
	 			delete [] pmatrix;
			return EXIT_FAILURE;
		}catch(...){
			std::cout << "Unknown failure in trying to load covariates matrix.\n";
	 		delete [] aux;
	 		delete [] evectors;
	 		delete [] h_y;
	 		if(get_pval)
	 			delete [] pmatrix;
			return EXIT_FAILURE;
		}
		if(cov == NULL){
			delete [] aux;
			delete [] evectors;
			delete [] h_y;
	 		if(get_pval)
	 			delete [] pmatrix;
			return EXIT_FAILURE;
		}


		if(n_subjects_1 != n_subjects){
			printf("Rows of trait matrix are not equal to rows of covariate matrix.\n  The number of rows in the trait matrix is %i and the number of rows in the covariate matrix is %i.\n",
					n_subjects, n_subjects_1);
			delete [] aux;
			delete [] evectors;
			delete [] h_y;
			delete [] cov;
	 		if(get_pval)
	 			delete [] pmatrix;
			return EXIT_FAILURE;
		}



 	}



	size_t qsub_batch_size;


	if(use_qsub){
		qsub_batch_size = n_voxels/total_nodes;

		if(node_number + 1 == total_nodes  && n_voxels % total_nodes != 0){
			qsub_batch_size = n_voxels % total_nodes;
		}

	}else{
		qsub_batch_size = n_voxels;
	}




 	float * indicator;
 	float * pvals;
 	float * h2;
 	float  * conn_h2;
 	float * conn_indicator;
 	float * conn_pvalues;


 	if(get_connectivity){
 		conn_h2 = new float[qsub_batch_size*n_voxels];
 		conn_indicator = new float[qsub_batch_size*n_voxels];
 		if(get_pval){
 			conn_pvalues = new float[n_voxels*qsub_batch_size];
 		}
 	}else{
 	 	indicator = new float[n_voxels];
 	 	pvals = new float[n_voxels];
 	 	h2 = new float[n_voxels];
 	}

 	int error;
 	if(use_qsub){



 		error = call_cudafphi_qsub( conn_h2, conn_indicator,
 				  conn_pvalues,  h_y, aux,  cov,
 		            evectors, pmatrix,  covariates,  get_pval,  get_connectivity,  node_number,
 		             n_voxels,  n_subjects,  n_permutations,  qsub_batch_size,  (n_voxels/total_nodes),
 		             n_covariates,  selected_devices);



 		error = write_to_hdf5_qsub(file_out, conn_h2, conn_indicator , conn_pvalues ,qsub_batch_size, n_voxels/total_nodes, node_number,  n_voxels,  get_pval);





 	}else{

 		error = call_cudafphi(h2, indicator, pvals, conn_h2, conn_indicator, conn_pvalues, h_y, aux, cov, evectors, pmatrix, covariates, get_pval, get_connectivity, n_voxels, n_subjects, n_permutations, n_covariates, selected_devices);
 	 	if(get_connectivity){
 	 		error = write_to_hdf5_qsub(file_out, conn_h2, conn_indicator , conn_pvalues ,qsub_batch_size, n_voxels, n_voxels,  n_voxels,  get_pval);
 		}else{
 	 	 	write_to_file(file_out,header_filename, h2, indicator, pvals, get_pval, n_voxels);
 		}
 	}
 	delete [] aux;
 	delete [] evectors;
 	delete [] h_y;
 	if(covariates)
 		delete [] cov;
 	if(get_pval)
 		delete [] pmatrix;




 	if(get_connectivity == false){


 	 	delete [] h2;
 	 	if(get_pval)
 	 		delete [] pvals;
 	 	delete [] indicator;



 	}else{
 		delete [] conn_h2;
 		delete [] conn_indicator;
 		if(get_pval){
 			delete [] conn_pvalues;
 		}
 	}


 	if(error == 0)
 		return EXIT_FAILURE;

 	if(debug)
 		close(error_file);
 	return EXIT_SUCCESS;
}




