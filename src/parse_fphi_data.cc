#include <unordered_map>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <algorithm>
#include <string>
#include <iterator>
#include "solar.h"
#include <fstream>
#include <list>
using namespace std;




static unordered_map<string, list<float>  > trait_map;
static unordered_map<string, list<int8_t>  > snp_map;
static int n_snps;
static unordered_map<string, list<float> > covariate_map;
static vector<int> ibdids;
static list<string> phenotype_ids;
static list<string> ids;
static list<int> trait_indices_list;
static list<int> cov_indices_list;

extern "C" void symeig_ (int*, double*, double*, double*, double*, int*);
bool is_bad_char(char c){
	if(isgraph(c) || c == '\n'){
		return false;
	}else{
		return true;
	}
}


void calculate_eigenvectors (double * phi2, double * eigenvalues, double * eigenvectors ,int n)
{
    double* e =  new double[n];
    memset(e, 0, sizeof(double)*n);
    int * info = new int;
    *info  = 0;
    symeig_(&n, phi2, eigenvalues, e, eigenvectors, info);
    delete [] e;
    delete [] info;
}

static void append_ibdids_and_ids(string id, int ibdid){

	list<string>::iterator iter = find(phenotype_ids.begin(), phenotype_ids.end(), id);
	if(iter != phenotype_ids.end()){
		ibdids.push_back(ibdid);
		ids.push_back(id);
		phenotype_ids.erase(iter);
	}
}


/*
static void remove_row(Eigen::MatrixXd& matrix,  int rowToRemove)
{
     int numRows = matrix.rows()-1;
     int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void remove_column(Eigen::MatrixXd& matrix,  int colToRemove)
{
     int numRows = matrix.rows();
     int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

int create_eigenvectors_and_auxiliary_matrix_gwas(Tcl_Interp * interp, Eigen::MatrixXd & Hx_aux,Eigen::MatrixXd & aux, Eigen::MatrixXd & Hx_eigenvectors, Eigen::MatrixXd & eigenvectors, \ 
		Eigen::MatrixXd & trait_matrix,Eigen::MatrixXd & snp_matrix, Eigen::MatrixXd & cov_matrix){
	
	Matrix* pm2;

	pm2 = Matrix::find("phi2");
	if (!pm2) {
	    Solar_Eval(interp, "loadkin");
	    pm2 = Matrix::find("phi2");
	    if(!pm2){
			RESULT_LIT("phi2 matrix could not be loaded.\n");
			return TCL_ERROR;
		}
	}	

/*		Eigen::MatrixXd phi2_matrix(ibdids.size(), ibdids.size());
	double phi2_value;
	
		for(int col_idx = 0; col_idx < ibdids.size(); col_idx++){
			for(int row_idx = col_idx; row_idx < ibdids.size(); row_idx++) {
				phi2_value = pm2->get(ibdids[row_idx], ibdids[col_idx]);
				phi2_matrix(row_idx, col_idx) = phi2_value;
				phi2_matrix(col_idx, row_idx) = phi2_value;
		//		phi2_matrix[col_idx*ibdids.size() + row_idx] = phi2_value;
		//		phi2_matrix[row_idx*ibdids.size() + col_idx] = phi2_value;
			
			}
		
		}	

			
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(phi2_matrix);
	eigenvectors = es.eigenvectors().householderQr().householderQ();
	aux.col(0) = Eigen::ArrayXd::Ones(ibdids.size()).matrix();
	aux.col(1) = es.eigenvalues();	
	
	
//	double*  phi2_matrix = new double[ibdids.size()*ibdids.size()];
	Eigen::MatrixXd phi2_matrix(ibdids.size(), ibdids.size());
	try{
		for(int col_idx = 0; col_idx < ibdids.size(); col_idx++){
			for(int row_idx = col_idx; row_idx < ibdids.size(); row_idx++) {
			
				double phi2_value = pm2->get(ibdids[row_idx], ibdids[col_idx]);
				phi2_matrix(row_idx, col_idx) = phi2_value;
				phi2_matrix(col_idx, row_idx) = phi2_value;
		//		phi2_matrix[col_idx*ibdids.size() + row_idx] = phi2_value;
		//		phi2_matrix[row_idx*ibdids.size() + col_idx] = phi2_value;
			
			}
		}	
	} catch(...) {
		RESULT_LIT("Error in reading in the phi2 matrix.\n");
		return TCL_ERROR;
	}	
	

	Hx_aux.resize(phi2_matrix.rows(), 2);	

	if(cov_matrix.cols() != 0){
		Eigen::MatrixXd hat_x  = (Eigen::MatrixXd::Identity(cov_matrix.rows(), cov_matrix.rows())  - \
							cov_matrix*(cov_matrix.transpose()*cov_matrix).inverse()\
							*cov_matrix.transpose());
		Eigen::MatrixXd Hx_phi2_matrix = hat_x*phi2_matrix*hat_x;
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_1(Hx_phi2_matrix);
		

		Hx_aux.col(0) = Eigen::ArrayXd::Ones(ibdids.size()).matrix();
		Hx_aux.col(1) = es_1.eigenvalues();	
		Hx_eigenvectors = es_1.eigenvectors();		
		double evalue;
		
		for(int idx = 0 ; idx < Hx_aux.rows(); idx++){
			evalue = Hx_aux(idx, 1);				
			if(evalue == 0.0){
				remove_row(Hx_aux, idx);
				remove_row(trait_matrix, idx);
				remove_row(Hx_eigenvectors, idx);
				remove_column(Hx_eigenvectors, idx);
				remove_row(phi2_matrix, idx);
				remove_column(phi2_matrix, idx);				
				remove_row(snp_matrix, idx);
				remove_row(cov_matrix, idx);
				idx--;
			}		
		}
		Hx_eigenvectors = Hx_eigenvectors.householderQr().householderQ();
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_2(phi2_matrix);
		aux.resize(phi2_matrix.rows(), 2);	
		aux.col(0) = Eigen::ArrayXd::Ones(ibdids.size()).matrix();
		aux.col(1) = es_2.eigenvalues();
		eigenvectors = es_2.eigenvectors().householderQr().householderQ();		
									
	}else{
		
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(phi2_matrix);
		aux.resize(phi2_matrix.rows(), 2);	
		aux.col(0) = Eigen::ArrayXd::Ones(ibdids.size()).matrix();
		aux.col(1) = es.eigenvalues();
		eigenvectors = es.eigenvectors().householderQr().householderQ();	
		Hx_eigenvectors = eigenvectors;
		Hx_aux = aux;		
		
	}
	
	return TCL_OK;
}*/



int create_eigenvectors_and_auxiliary_matrix(Tcl_Interp * interp, Eigen::MatrixXd & aux, Eigen::MatrixXd & eigenvectors){
	
	Matrix* pm2;

	pm2 = Matrix::find("phi2");
	if (!pm2) {
	    Solar_Eval(interp, "loadkin");
	    pm2 = Matrix::find("phi2");
	    if(!pm2){
			cout << "phi2 matrix could not be loaded.\n";
			return TCL_ERROR;
		}
	}	
	cerr << "Performing EVD...\n";
	double * phi2_matrix = new double[ids.size()*ids.size()];
	int n_subjects = ids.size();
	double phi2_value;
	
		for(int col_idx = 0; col_idx < n_subjects; col_idx++){
			for(int row_idx = col_idx; row_idx < n_subjects; row_idx++) {
				phi2_value = pm2->get(ibdids[row_idx], ibdids[col_idx]);
				phi2_matrix[col_idx*n_subjects + row_idx] = phi2_value;
				phi2_matrix[row_idx*n_subjects + col_idx] = phi2_value;
			
			}
		
		}	
	double * d_eigenvectors = new double[ids.size()*ids.size()];
	double * d_eigenvalues = new double[ids.size()];

	calculate_eigenvectors (phi2_matrix, d_eigenvalues, d_eigenvectors , ids.size());
	
	delete [] phi2_matrix;
	
	aux.resize(ids.size(),2);
	eigenvectors.resize(ids.size(), ids.size());

	for(int row = 0; row < n_subjects; row++){
		aux(row, 0) = 1.f;
		aux(row, 1) = d_eigenvalues[row];
		for(int col = 0 ; col < n_subjects; col++){
			eigenvectors(row, col)= d_eigenvectors[col*n_subjects + row];
		}
	}
	delete [] d_eigenvectors;
	delete [] d_eigenvalues;
	
	return TCL_OK;
}	
static int create_snps_matrix(Tcl_Interp * interp, Eigen::MatrixXd & snps_matrix){
	list<int8_t> snp_list;
	list<int8_t>::iterator snp_iter;
	int index;
	int row_index = 0;
	list<double> snp_means;
	list<double>::iterator mean_iter;
	list<int> snp_count;
	list<int>::iterator count_iter;
	int8_t current_snp;
	for(list<string>::iterator id_iter = ids.begin(); id_iter != ids.end(); id_iter++, row_index++){
		
		snp_list = snp_map.at(*id_iter);
		snp_iter = snp_list.begin();
		if(row_index == 0){
			snps_matrix.resize(ids.size(), snp_list.size());
			snp_means.resize(snp_list.size());
			snp_count.resize(snp_list.size());
			count_iter = snp_count.begin();
			for(mean_iter = snp_means.begin(); mean_iter != snp_means.end(); mean_iter++, count_iter++){
				*mean_iter = 0.0;
			}
		}
		index = 0;
		mean_iter = snp_means.begin();
		count_iter = snp_count.begin();
		for(snp_iter = snp_list.begin(); snp_iter != snp_list.end(); snp_iter++, index++, mean_iter++, count_iter++){
			if((current_snp = *snp_iter) != -1){
				*mean_iter += current_snp/2.0;
				*count_iter = *count_iter + 1;
				snps_matrix(row_index, index) = current_snp/2.0;
			}else{
				snps_matrix(row_index, index) = -1.0;
			}
			
		}
		
		
	}
	mean_iter = snp_means.begin();
	count_iter = snp_count.begin();
	for(int col_index = 0; col_index < snps_matrix.cols() ;col_index++, mean_iter++, count_iter++){
		*mean_iter /= *count_iter;
		for(int row_index = 0; row_index< snps_matrix.rows(); row_index++){
			if(snps_matrix(row_index,col_index) != -1.0){
				snps_matrix(row_index, col_index) -= *mean_iter; 
			}else{
				snps_matrix(row_index,col_index) = 0.0;
			}
		}
	}
		
	
	
	return TCL_OK;
	
}


static int create_trait_and_covariate_matrix(Tcl_Interp * interp, Eigen::MatrixXd & trait_matrix, Eigen::MatrixXd & covariate_matrix, \
									vector<string> unique_cov_list, vector<string> cov_list, int n_traits){									
	list<float>  traits;
	list<float> covariates;
	int row_idx = 0;
	string id;
	list<float>::iterator row_iter;
	Eigen::MatrixXd unique_cov_matrix;
	if(unique_cov_list.size() != 0)  unique_cov_matrix.resize(ids.size(), unique_cov_list.size());
	for(list<string>::iterator id_iter = ids.begin() ; id_iter != ids.end(); id_iter++){
		id = *id_iter;
		traits = trait_map.at(id);
		if(cov_list.size() != 0) covariates = covariate_map.at(id);
		row_iter = traits.begin();
		for( int trait = 0 ; trait < n_traits; trait++, row_iter++)  \
		trait_matrix(row_idx, trait) = *row_iter;
		
		row_iter = covariates.begin();
		for( int cov = 0 ; cov < unique_cov_list.size(); cov++, row_iter++) \
			unique_cov_matrix(row_idx, cov) = *row_iter;
			
		row_idx++;
	}
		
	for(int col = 0; col < unique_cov_matrix.cols(); col++){
	
		unique_cov_matrix.col(col) =unique_cov_matrix.col(col).array() - unique_cov_matrix.col(col).mean();
	}
	
	for(int cov = 0; cov < cov_list.size(); cov++){
		vector<string> cov_terms;
		vector<int> cov_indices;
		string cov_string = string(cov_list[cov]);
		string cov_term;
		char c;
		for(int i = 0 ; i < cov_string.size() ; i++){
			c = cov_string[i];
			if(c == '*'){
				cov_terms.push_back(cov_term);
				cov_term = "";
			}else{
				cov_term += c;
			}
		}
		cov_terms.push_back(cov_term);
		
		for(vector<string>::iterator cov_iter = cov_terms.begin(); cov_iter != cov_terms.end(); cov_iter++){
			int term_idx = 0;
			
			while(cov_iter->compare(unique_cov_list[term_idx]) != 0) term_idx++;
			
			if(cov_iter == cov_terms.begin())
				covariate_matrix.col(cov) = unique_cov_matrix.col(term_idx);
			else
				covariate_matrix.col(cov) = (unique_cov_matrix.col(term_idx).array()*covariate_matrix.col(cov).array()).matrix();
		}
	}
	
	if(covariate_matrix.cols() != 0){
			
		Eigen::MatrixXd cov_matrix_with_intercept(covariate_matrix.rows(), covariate_matrix.cols() + 1);
		
		for(int col = 0; col < covariate_matrix.cols() + 1; col++){
			if(col != covariate_matrix.cols()){
				cov_matrix_with_intercept.col(col) = (covariate_matrix.col(col).array() - covariate_matrix.col(col).mean())/(covariate_matrix.col(col).norm()/sqrt(covariate_matrix.rows() - 1));
			}else{
				cov_matrix_with_intercept.col(col) = Eigen::ArrayXd::Ones(covariate_matrix.rows()).matrix();
			}
		}	
		covariate_matrix = cov_matrix_with_intercept;
	}
		
	/*if(covariate_matrix.cols() != 0){
			
		Eigen::MatrixXd cov_matrix_with_intercept(covariate_matrix.rows(), covariate_matrix.cols());
		
		for(int col = 0; col < covariate_matrix.cols(); col++){
			
				cov_matrix_with_intercept.col(col) = (covariate_matrix.col(col).array() - covariate_matrix.col(col).mean())/(covariate_matrix.col(col).norm()/sqrt(covariate_matrix.rows() - 1));
			
		}	
		covariate_matrix = cov_matrix_with_intercept;
	}
*/
	return TCL_OK;
}


	

static bool is_whitespace(char c){
	if(c == '\t' || c == ' '){
		return true;
	}
	
	return false;
}
			
static bool not_is_graph( const char  c) { return !isgraph(c); }

static vector<string> parse_line_vector(  string   line, bool to_upper){
	
	if(to_upper)
		transform(line.begin(), line.end(), line.begin(), ::toupper);
	vector<string> parsed_line;
	remove_if(line.begin() ,line.end(), not_is_graph);
	string token;	
   

    for(string::iterator iter = line.begin(); iter != line.end(); iter++){
		if(*iter != ','){
			token += *iter;
		}else{
			parsed_line.push_back(token);
			token.clear();
		}
		
	} 
    parsed_line.push_back(token);
	return parsed_line;
}

static vector<string> parse_line_to_max_index(  string   line, int max_index){
	

	vector<string> parsed_line;
	remove_if(line.begin() ,line.end(), not_is_graph);
	string token;	
   
	int index = 0;
    for(string::iterator iter = line.begin(); iter != line.end(); iter++){
		if(*iter != ','){
			token += *iter;
		}else{
			parsed_line.push_back(token);
			token.clear();
			if(++index > max_index)
				break;
		}
		
	} 
	if(index <= max_index)
		parsed_line.push_back(token);
	return parsed_line;
}

bool is_nan(double val) { if(val == val) return false; return true; }
static list<string> parse_line_list(string  line, bool to_upper){
	
	if(to_upper)
		transform(line.begin(), line.end(), line.begin(), ::toupper);
	remove_if(line.begin() ,line.end(), not_is_graph);
	string token;	
    list<string> parsed_line;
    
    for(string::iterator iter = line.begin(); iter != line.end(); iter++){
		if(*iter != ','){
			token += *iter;
		}else{
			parsed_line.push_back(token);
			token.clear();
		}
		
	} 
    parsed_line.push_back(token);

	
	return parsed_line;
}	
static int create_fphi_id_list(Tcl_Interp * interp, string phenotype_filename, vector<string> traits, vector<string> covs){
	
	cerr << "Reading phenotype file...\n";
	
	const char * errmsg = 0;
	SolarFile * file = SolarFile::open("fphi", phenotype_filename.c_str(), &errmsg);
	
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	int field_count;
	const char ** fields = file->names(&field_count, &errmsg);
	
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	
	file->start_setup(&errmsg);
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	file->setup("id", &errmsg);
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}	
	char ** file_data;
	while (0 != (file_data = file->get (&errmsg))){
		
		phenotype_ids.push_back(string(file_data[0]));
		
		
	}
	file->rewind(&errmsg);	
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	file->start_setup(&errmsg);
	
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}

	
	for(vector<string>::iterator trait_iter = traits.begin(); trait_iter != traits.end(); trait_iter++){
		file->setup(trait_iter->c_str(), &errmsg);
		if(errmsg){
			RESULT_LIT(errmsg);
			return TCL_ERROR;
		}		
	}
	
	
	list<float> row_values(traits.size());
	list<float>::iterator row_iter;
	size_t index;
	size_t n_traits = traits.size();
	bool skip_row;
	list<string>::iterator id_iter = phenotype_ids.begin();
	list<string> bad_ids;
	while (0 != (file_data = file->get (&errmsg))){
		
		row_iter = row_values.begin();
		skip_row = false;
		for(index = 0; index < n_traits; index++, row_iter++){
			if(StringCmp(file_data[index], 0, case_ins)){
				*row_iter = strtof(file_data[index], NULL);
			}else{
				skip_row = true;
				break;
			}
		}
		
		if(skip_row){
			bad_ids.push_back(*id_iter++);
			continue;
		}
		trait_map[*id_iter++] = row_values;
	}
	file->rewind(&errmsg);	
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	file->start_setup(&errmsg);
	
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	
	
	for(vector<string>::iterator cov_iter = covs.begin(); cov_iter != covs.end(); cov_iter++){
		file->setup(cov_iter->c_str(), &errmsg);
		if(errmsg){
			RESULT_LIT(errmsg);
			return TCL_ERROR;
		}		
	}
	int n_covs = covs.size();
	row_values.resize(n_covs);
	
	id_iter = phenotype_ids.begin();
	while (0 != (file_data = file->get (&errmsg))){
		
		row_iter = row_values.begin();
		skip_row = false;
		for(index = 0; index < n_covs; index++, row_iter++){
			if(!StringCmp(file_data[index], "F", case_ins)){
				*row_iter = 2;
			
			}else if(!StringCmp(file_data[index], "M", case_ins)){
				*row_iter = 1;
			}else if(StringCmp(file_data[index], 0, case_ins)){
				*row_iter = strtof(file_data[index], NULL);
			}else{
				skip_row = true;
				break;
			}
		}
		
		if(skip_row){
			id_iter = phenotype_ids.erase(id_iter);
			continue;
		}
		covariate_map[*id_iter++] = row_values;
	}
	
	for(list<string>::iterator bad_id_iter = bad_ids.begin(); bad_id_iter != bad_ids.end(); bad_id_iter++){
		id_iter = find(phenotype_ids.begin(), phenotype_ids.end(), *bad_id_iter);
		if(id_iter != phenotype_ids.end()){
			phenotype_ids.erase(id_iter);
		}
	}
	
	SolarFile * ped_file = SolarFile::open("fphi", "pedindex.out", &errmsg);
		if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	ped_file->start_setup(&errmsg);
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	ped_file->setup("id", &errmsg);
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
    int ibdid = 1;
	while (0 != (file_data = ped_file->get (&errmsg))){	
		append_ibdids_and_ids(string(file_data[0]), ibdid++);		
	}				
				
	return TCL_OK;
}	
static int create_gwas_id_list(Tcl_Interp * interp, string phenotype_filename, vector<string> traits, vector<string> covs, vector<string> & snp_names,
				int n_traits, int n_covs){
	cerr << "Reading phenotype file...\n";				
	const char * errmsg = 0;
	SolarFile * file = SolarFile::open("gwas", phenotype_filename.c_str(), &errmsg);
	
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	int field_count;
	const char ** fields = file->names(&field_count, &errmsg);
	
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	
	for(int field = 0; field < field_count ; field++){
		//cout << fields[field] << endl;
		if(strstr(Strdup(fields[field]), "snp_")){
			snp_names.push_back(string(fields[field]));
		}
	}	
	file->start_setup(&errmsg);
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	file->setup("id", &errmsg);
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}	
	char ** file_data;
	
	while (0 != (file_data = file->get (&errmsg))){
		
		phenotype_ids.push_back(string(file_data[0]));
		
		
	}
	file->rewind(&errmsg);	
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	file->start_setup(&errmsg);
	
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
		
	for(vector<string>::iterator trait = traits.begin(); trait != traits.end() ; trait++){
		file->setup(trait->c_str(), &errmsg);				
		if(errmsg){
			RESULT_LIT(errmsg);
			return TCL_ERROR;
		}
	}
	list<float> row_values(traits.size());
	list<string>::iterator id_iter = phenotype_ids.begin();
	list<float>::iterator row_iter;
	bool skip_row;
	list<string> bad_ids;
	
	while (0 != (file_data = file->get (&errmsg))){
		row_iter = row_values.begin();
		skip_row = false;
		for(int idx = 0; idx < traits.size(); idx++, row_iter++){
			if(StringCmp(file_data[idx],0, 0)){
				*row_iter = strtof((const char*)file_data[idx], NULL);
			}else{
				skip_row = true;
				break;
			}
		}
		if(skip_row){
			bad_ids.push_back(*id_iter);
			id_iter++;	
			continue;
		}
		
		trait_map[*id_iter] = row_values;
		id_iter++;
	}	
	string str;

	file->rewind(&errmsg);	
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	file->start_setup(&errmsg);
	
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
						
	for(vector<string>::iterator snp_iter = snp_names.begin(); snp_iter != snp_names.end(); snp_iter++){
		file->setup(snp_iter->c_str(), &errmsg);
	}
	id_iter = phenotype_ids.begin();
//	row_values.resize(snp_names.size());
	list<string>::iterator bad_id_iter;
	
	list<int8_t> snp_row_values(snp_names.size());
	list<int8_t>::iterator snp_iter;
	while (0 != (file_data = file->get (&errmsg))){	
		snp_iter = snp_row_values.begin();
		skip_row = false;
		for(int idx = 0; idx < snp_names.size(); idx++, snp_iter++){
			if(StringCmp(file_data[idx],0, 0)){
				*snp_iter = strtol((const char*)file_data[idx], NULL, 10);
			}else{
				*snp_iter = -1;
			}
		}
		
		if ( all_of(snp_row_values.begin(), snp_row_values.end(), [](int8_t i){return i == -1;}) ){
			id_iter = phenotype_ids.erase(id_iter);
			continue;
		}		
		snp_map[*id_iter] = snp_row_values;
		id_iter++;
	}
	for(bad_id_iter = bad_ids.begin(); bad_id_iter != bad_ids.end(); bad_id_iter++){
		id_iter = find(phenotype_ids.begin(), phenotype_ids.end(), *bad_id_iter);
		if(id_iter != phenotype_ids.end())
			phenotype_ids.erase(id_iter);
	}
	
	SolarFile * ped_file = SolarFile::open("gwas", "pedindex.out", &errmsg);
		if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	ped_file->start_setup(&errmsg);
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	ped_file->setup("id", &errmsg);
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
    int ibdid = 1;
	while (0 != (file_data = ped_file->get (&errmsg))){	
		append_ibdids_and_ids(string(file_data[0]), ibdid++);		
	}	
	
	return TCL_OK;
}


int generate_gwas_matrices(Tcl_Interp * interp, Eigen::MatrixXd  &trait_vector, Eigen::MatrixXd  & cov_matrix,Eigen::MatrixXd & l_Hx_aux, Eigen::MatrixXd & l_aux, Eigen::MatrixXd & l_Hx_eigenvectors,\
Eigen::MatrixXd & l_eigenvectors, Eigen::MatrixXd & snps_matrix, vector<string> & snp_names, string * trait_names, int n_traits){
					
	static Eigen::MatrixXd eigenvectors;
	static Eigen::MatrixXd aux;
	static Eigen::MatrixXd Hx_aux;
	static Eigen::MatrixXd Hx_eigenvectors;
	static vector<string> saved_cov_list;
	static list<string> saved_ids;					

	vector<string> cov_list;		
	vector<string> unique_cov_terms;	
	int success;				

	vector<string> string_vector(n_traits);
	string current_trait;
	for(int i = 0; i <n_traits; i++){
		current_trait = trait_names[i];
		transform(current_trait.begin(), current_trait.end(), current_trait.begin(), ::toupper);
		string_vector[i] = current_trait;
	}
	Covariate * c;
	for (int i = 0;( c = Covariate::index(i)); i++)
	{
		char buf[512];
		size_t first_pos = 0;

		c->fullname (&buf[0]);
		string current_cov = string(&buf[0]);			
		transform(current_cov.begin(), current_cov.end(), current_cov.begin(), ::toupper);
		cov_list.push_back(current_cov.c_str());	
	}
	for(int cov = 0 ; cov < cov_list.size() ; cov++){
		string current_cov = cov_list[cov];
		string current_term;
		for(string::iterator it = current_cov.begin(); it != current_cov.end(); it++){
			if(*it == '*'){
				if(find(unique_cov_terms.begin(), unique_cov_terms.end(), current_term) == unique_cov_terms.end()) \
				unique_cov_terms.push_back(current_term);
				current_term = "";
			}else{
				current_term += *it;
			}
		}
		if(find(unique_cov_terms.begin(), unique_cov_terms.end(), current_term) == unique_cov_terms.end()) unique_cov_terms.push_back(current_term);
	}
	

	string phenotype_filename = Phenotypes::filenames();
	vector<string> cov_terms_copy(unique_cov_terms);
	success = create_gwas_id_list(interp,  phenotype_filename, string_vector, cov_terms_copy, snp_names, n_traits, unique_cov_terms.size());
	
	if(success == TCL_ERROR){
		trait_map.clear();
		covariate_map.clear();
		ibdids.clear();
		phenotype_ids.clear();
		ids.clear();
		return success;
	 }
 
	trait_vector.resize(ids.size(), n_traits);
		

	cov_matrix.resize(ids.size(), cov_list.size());
	
	success = create_trait_and_covariate_matrix(interp,  trait_vector, cov_matrix, \
									unique_cov_terms, cov_list,  n_traits);		
									

								
	if(success == TCL_ERROR){
		trait_map.clear();
		covariate_map.clear();
		ibdids.clear();
		phenotype_ids.clear();
		ids.clear();
		snp_map.clear();
		return success;
	 }	
	 	
	
	success = create_snps_matrix(interp, snps_matrix);
	if(success == TCL_ERROR){
		trait_map.clear();
		covariate_map.clear();
		ibdids.clear();
		phenotype_ids.clear();
		ids.clear();
		snp_map.clear();
		return success;
	 }		
	if(saved_ids != ids || saved_cov_list != cov_list){
		
		success = create_eigenvectors_and_auxiliary_matrix(interp, aux ,eigenvectors);
		
		if(success == TCL_ERROR){
			trait_map.clear();
			covariate_map.clear();
			ibdids.clear();
			phenotype_ids.clear();
			ids.clear();
			snp_map.clear();
			return success;
		}
	}
		
		
	saved_ids = ids;
	saved_cov_list = cov_list;
	trait_map.clear();
	covariate_map.clear();
	ibdids.clear();
	phenotype_ids.clear();
	ids.clear();
	snp_map.clear();
	l_aux = aux;
	l_eigenvectors = eigenvectors;
	l_Hx_aux = Hx_aux;
	l_Hx_eigenvectors = Hx_eigenvectors;
	return TCL_OK;
}	


int generate_fphi_matrices(Tcl_Interp * interp, Eigen::MatrixXd  & l_trait_vector, Eigen::MatrixXd  & l_cov_matrix, Eigen::MatrixXd & l_aux, Eigen::MatrixXd & l_eigenvectors,\
				vector<string> trait_names, int n_traits){
					
	static Eigen::MatrixXd eigenvectors;
	static Eigen::MatrixXd aux;
	static list<string> saved_ids;					
	Eigen::MatrixXd trait_vector;
	Eigen::MatrixXd cov_matrix;
							
	vector<string> cov_list;		
	vector<string> unique_cov_terms;	
	int success;				

	Covariate * c;
	for (int i = 0;( c = Covariate::index(i)); i++)
	{
		char buf[512];
		size_t first_pos = 0;

		c->fullname (&buf[0]);
		string current_cov = string(&buf[0]);			
		transform(current_cov.begin(), current_cov.end(), current_cov.begin(), ::toupper);
		cov_list.push_back(current_cov.c_str());	
	}
	for(int cov = 0 ; cov < cov_list.size() ; cov++){
		string current_cov = cov_list[cov];
		string current_term;
		for(string::iterator it = current_cov.begin(); it != current_cov.end(); it++){
			if(*it == '*'){
				if(find(unique_cov_terms.begin(), unique_cov_terms.end(), current_term) == unique_cov_terms.end()) \
				unique_cov_terms.push_back(current_term);
				current_term = "";
			}else{
				current_term += *it;
			}
		}
		if(find(unique_cov_terms.begin(), unique_cov_terms.end(), current_term) == unique_cov_terms.end()) unique_cov_terms.push_back(current_term);
	}
	
	
	//success = create_id_list(new_trait, unique_cov_terms, new_ibdids, new_ids);

	//if(success == TCL_ERROR) return success;

	string phenotype_filename = Phenotypes::filenames();
	vector<string> cov_terms_copy(unique_cov_terms);
	success = create_fphi_id_list(interp,  phenotype_filename, trait_names, unique_cov_terms);
	
	if(success == TCL_ERROR){
		trait_map.clear();
		covariate_map.clear();
		ibdids.clear();
		phenotype_ids.clear();
		ids.clear();
		return success;
	 }
	trait_vector.resize(ids.size(), n_traits);
		

	cov_matrix.resize(ids.size(), cov_list.size());
	

	success = create_trait_and_covariate_matrix(interp,  trait_vector, cov_matrix, \
									unique_cov_terms, cov_list,  n_traits);		
									

								
	if(success == TCL_ERROR){
		trait_map.clear();
		covariate_map.clear();
		ibdids.clear();
		phenotype_ids.clear();
		ids.clear();
		return success;
	 }	
	
	if(saved_ids != ids){
		eigenvectors.resize(ibdids.size(), ibdids.size());
		aux.resize(ibdids.size(), 2);	
		success = create_eigenvectors_and_auxiliary_matrix(interp, aux, eigenvectors);
		if(success == TCL_ERROR){
			trait_map.clear();
			covariate_map.clear();
			ibdids.clear();
			phenotype_ids.clear();
			ids.clear();
			return success;
		}
	}
		
		
	saved_ids = ids;
	trait_map.clear();
	covariate_map.clear();
	ibdids.clear();
	phenotype_ids.clear();
	ids.clear();
	l_trait_vector = trait_vector;	
	l_eigenvectors = eigenvectors;
	l_cov_matrix = cov_matrix;
	l_aux = aux;
	

	return TCL_OK;
}


