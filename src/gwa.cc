#include "solar.h"
#define EIGEN_MAX_CPP_VER 11
#include "Eigen/Dense"
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include "plinkio.h"
#include <unordered_map>
#include <list>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <omp.h>
using namespace std;

static const int GWAS_BATCH_SIZE = 10000;

void calculate_eigenvectors (double * phi2, double * eigenvalues, double * eigenvectors ,int n);
void calculate_sigma_A_and_sigma_E_fphi(Eigen::MatrixXd Y, Eigen::MatrixXd Z, Eigen::MatrixXd ZTZI, Eigen::VectorXd & Sigma);
//int generate_gwas_matrices(Tcl_Interp * interp, Eigen::MatrixXd  &trait_vector, Eigen::MatrixXd  & cov_matrix,Eigen::MatrixXd & l_Hx_aux, Eigen::MatrixXd & l_aux, Eigen::MatrixXd & l_Hx_eigenvectors,\
Eigen::MatrixXd & l_eigenvectors, Eigen::MatrixXd & snps_matrix,vector<string> & snp_list, string * trait_names, int n_traits);
static int generate_gwas_matrices(Tcl_Interp * interp, Eigen::MatrixXd & trait_vector, int8_t ** snp_values, Eigen::MatrixXd & evectors, Eigen::MatrixXd & aux, 
		string trait_name, vector<string> & snp_list, string plink_filename, string phenotype_filename, int & n_subjects){
	static list<string> saved_ids;
	static Eigen::MatrixXd saved_evectors;
	static Eigen::MatrixXd saved_aux;		
	const char * errmsg = 0;
	cerr << "Reading phenotype file...\n";
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
	list<string> phenotype_ids;

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
	file->setup(trait_name.c_str(), &errmsg);
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	unordered_map<string , double > trait_map;
	double row_value;
	list<string>::iterator id_iter = phenotype_ids.begin();
	bool skip_row;
	list<string> bad_ids;

	while (0 != (file_data = file->get (&errmsg))){
		
		if(StringCmp(file_data[0], 0 , 0)){
			trait_map[*id_iter] = strtod((const char*)file_data[0], NULL);
			id_iter++;
		}else{
			bad_ids.push_back(*id_iter++);
		}
	}	
	string str;
	pio_file_t plink_file;
	int sample_size;
	unsigned int n_snps;
	unordered_map<string, int> plink_map;	
	list<int8_t> snp_row_values;
	unordered_map<string, list<int8_t> > snp_map;
	unsigned int snp_idx;
	list<string> plink_ids;	
	string snp_name;	
	if(plink_filename != ""){
		pio_status_t status = pio_open(&plink_file, plink_filename.c_str());
	

		if(status != PIO_OK){
			pio_close(&plink_file);
			if(status == P_FAM_IO_ERROR){
				RESULT_LIT("Error in loading .fam file");
				return TCL_ERROR;
			}else if (status == P_BIM_IO_ERROR){
				RESULT_LIT("Error in loading .bim file");
				return TCL_ERROR;
			}else if (status == P_BED_IO_ERROR){
				RESULT_LIT("Error in loading .bed file");
				return TCL_ERROR;
			}
		}
  
		n_snps  = plink_file.bed_file.header.num_loci;
		sample_size = plink_file.bed_file.header.num_samples;

	
	
		string id;
		pio_sample_t * sample;
		pio_locus_t * locus;
		snp_row_values.resize(n_snps);	
	
		for(snp_idx = 0; snp_idx < n_snps; snp_idx++){
			locus = bim_get_locus(&plink_file.bim_file, snp_idx);
			snp_name = locus->name;
			snp_list.push_back(snp_name);
		}
		for(list<string>::iterator bad_id_iter = bad_ids.begin(); bad_id_iter != bad_ids.end(); bad_id_iter++){
			id_iter = find(phenotype_ids.begin(), phenotype_ids.end(), *bad_id_iter);
			phenotype_ids.erase(id_iter);
		} 
		for(int subject = 0 ; subject < sample_size; subject++){
			
			sample = fam_get_sample(&plink_file.fam_file, subject); 
			id = sample->iid;
			id_iter = find(phenotype_ids.begin(), phenotype_ids.end(), id);
			if(id_iter == phenotype_ids.end()){
				id = sample->fid;
				id_iter = find(phenotype_ids.begin(), phenotype_ids.end(), id);
				if(id_iter == phenotype_ids.end()){
					continue;
				}
			}
			plink_map[id] = subject;
			plink_ids.push_back(id);
		

		}
		list<string>::iterator find_iter;
		for(id_iter = phenotype_ids.begin(); id_iter != phenotype_ids.end(); id_iter++){
			find_iter = find(plink_ids.begin(), plink_ids.end(), *id_iter);
			if(find_iter == phenotype_ids.end()){
				id_iter = phenotype_ids.erase(find_iter);
				id_iter--;
			}
		}		
	

		

	}else{
		
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
		for(unsigned int field = 0 ; field < field_count; field++){
			if(strstr(fields[field], "snp_")){
				snp_list.push_back(string(fields[field]));
			}
		}
		
		for(vector<string>::iterator snp_iter = snp_list.begin(); snp_iter != snp_list.end(); snp_iter++){
			file->setup(snp_iter->c_str(), &errmsg);
			*snp_iter = snp_iter->erase(0, 4);
		}
		snp_row_values.resize(snp_list.size());
		id_iter = phenotype_ids.begin();
		list<int8_t>::iterator snp_iter;
		while (0 != (file_data = file->get (&errmsg))){	
			snp_iter = snp_row_values.begin();
			skip_row = false;
			for(int idx = 0; idx < snp_list.size(); idx++, snp_iter++){
				if(StringCmp(file_data[idx],0, 0)){
					*snp_iter = strtol((const char*)file_data[idx], NULL, 10);
				}else{
					*snp_iter = 3;
				}
			}
			if ( all_of(snp_row_values.begin(), snp_row_values.end(), [](int8_t i){return i == 3;}) ){
				id_iter = phenotype_ids.erase(id_iter);
				continue;
			}
							
			snp_map[*id_iter] = snp_row_values;
			id_iter++;
		}
		
		for(list<string>::iterator bad_id_iter = bad_ids.begin(); bad_id_iter != bad_ids.end(); bad_id_iter++){
			id_iter = find(phenotype_ids.begin(), phenotype_ids.end(), *bad_id_iter);
			if(id_iter != phenotype_ids.end()){
				phenotype_ids.erase(id_iter);
			}
		} 					
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
	list<string> ids;
	vector<int> ibdids;
    int ibdid = 1;
	while (0 != (file_data = ped_file->get (&errmsg))){	
		id_iter = find(phenotype_ids.begin(), phenotype_ids.end(), string(file_data[0]));
		if(id_iter != phenotype_ids.end()){
			ids.push_back(*id_iter);
			ibdids.push_back(ibdid);
		}
		ibdid++;		
	}
	n_subjects = ids.size();	
	trait_vector.resize(n_subjects, 1);
	int trait_idx;
	int id_index = 0;
	for(id_iter = ids.begin(); id_iter != ids.end(); id_iter++, id_index++){
		row_value = trait_map[*id_iter];
		trait_vector(id_index, 0) = row_value;
	}
	trait_map.clear();
	if(saved_ids != ids){
		Matrix* pm2;

		pm2 = Matrix::find("phi2");
		if (!pm2) {
			Solar_Eval(interp, "loadkin");
			pm2 = Matrix::find("phi2");
			if(!pm2){
				RESULT_LIT("phi2 matrix could not be loaded");
				return TCL_ERROR;
			}
		}	

		double * phi2_matrix = new double[ids.size()*ids.size()];
		double phi2_value;
		cerr << "Performing EVD...\n";
	
			for(int col_idx = 0; col_idx < n_subjects; col_idx++){
				for(int row_idx = col_idx; row_idx < n_subjects; row_idx++) {
					phi2_value = pm2->get(ibdids[row_idx], ibdids[col_idx]);
					phi2_matrix[col_idx*n_subjects + row_idx] = phi2_value;
					phi2_matrix[row_idx*n_subjects + col_idx] = phi2_value;
			
				}
		
			}	
		double * d_eigenvectors = new double[ids.size()*ids.size()];
		double * eigenvalues = new double[ids.size()];
		calculate_eigenvectors (phi2_matrix, eigenvalues, d_eigenvectors , ids.size());
	
		delete [] phi2_matrix;
		aux.resize(n_subjects, 2);
		evectors.resize(n_subjects, n_subjects);
	
//	float * temp_evectors = new float[n_subjects*n_subjects];
//	float * temp_Z = new float[n_subjects*2];
	
		for(int row = 0; row < n_subjects; row++){
			aux(row, 0) = 1.0;
			aux(row, 1) = eigenvalues[row];
			for(int col = 0; col < n_subjects; col++){
				evectors(row, col) = d_eigenvectors[col*n_subjects + row];
			}
		}
		delete [] eigenvalues;
		delete [] d_eigenvectors;
		saved_ids = ids;
		saved_evectors = evectors;
		saved_aux = aux;
	}else{
		evectors = saved_evectors;
		aux = saved_aux;
	}
	n_snps = snp_list.size();	
	if(plink_filename != ""){
		snp_t * snp_buffer = new snp_t[sample_size];
		int8_t * temp_snps = new int8_t[n_snps*ids.size()];
		int8_t snp_value;
		pio_status_t status;
		cerr << "Reading plink file...\n";
		for(size_t snp = 0 ; snp < n_snps; snp++){
			status = pio_next_row(&plink_file, snp_buffer);
			if(status != PIO_OK){
				RESULT_LIT("Error in reading plink file");
				return TCL_ERROR;
			}		
			id_index = 0;
		
			for(list<string>::iterator id_iter = ids.begin(); id_iter != ids.end(); id_iter++, id_index++){
				temp_snps[snp*n_subjects + id_index] = snp_buffer[plink_map[*id_iter]];
			}
		
		
		}
		delete [] snp_buffer;
		*snp_values = temp_snps;
		pio_close(&plink_file);	
	}else{
		int8_t * temp_snps = new int8_t[n_snps*ids.size()];
		id_index = 0;
		list<int8_t>::iterator snp_iter;
		for(list<string>::iterator id_iter = ids.begin(); id_iter != ids.end(); id_iter++, id_index++){
			snp_row_values = snp_map[*id_iter];
			snp_iter = snp_row_values.begin();
			for(size_t snp = 0 ; snp < n_snps; snp++, snp_iter++){
				temp_snps[snp*n_subjects + id_index] = *snp_iter;
			}
		}
		*snp_values = temp_snps;
	}
	
	return TCL_OK;		
		
		
}				


static inline bool return_isspace (char  x){
	if(std::isspace(x))
		return true;
    else
	   return false;

}
extern int create_trait_vector(Tcl_Interp *interp);
static std::vector<string> parse_line(std::string currentLine, bool to_upper_case){

	currentLine.erase(std::remove_if(currentLine.begin(), currentLine.end(), return_isspace), currentLine.end());
	if(to_upper_case) transform(currentLine.begin(), currentLine.end(), currentLine.begin(), ::toupper);
	std::vector<string> output_vector;
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

static string parse_line_at_index(std::string currentLine, int index){

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

static int load_snps(vector<string> ids, Eigen::MatrixXd SNPs){

	string phenotype_file = string(Phenotypes::filenames());
	
	ifstream pheno_in(phenotype_file.c_str());
	
	string line;
	getline(pheno_in, line);
	
	vector<string> title = parse_line(line, true);
	
	int n_snps = 0;
	
	vector<string>::iterator id_iter =  find(title.begin(), title.end(), "ID");

		
	unsigned int id_index = distance(title.begin(), id_iter);	
	vector<int> snp_index;
	for(vector<string>::iterator iter = title.begin(); iter != title.end(); iter++) {
		string current_field = *iter;
		if(current_field.size() >= 4){
			if(current_field.substr(0, 4) == "SNP_"){
				n_snps++;
				snp_index.push_back(distance(title.begin(), iter));
			}
		}		
	}
	
	if(n_snps == 0){
		printf("No SNPs where detected in the phenotype file currently loaded.\nSNPs are detected by putting snp_ before the field name.\n");
		return TCL_ERROR;
	}
	string snp_str;
	SNPs.resize(ids.size(), snp_index.size());
	string id;
	vector<string>::iterator line_iter;
	int row_idx;
	while(getline(pheno_in, line)){
		

		
		bool continue_on = true;	
		vector<string> snp_vector;
		for(vector<int>::iterator snp_iter = snp_index.begin() ; snp_iter != snp_index.end(); snp_iter++){
			snp_str = parse_line_at_index(line, *snp_iter);
			if(snp_str == ""){
				continue_on =false;
				break;
			}
			snp_vector.push_back(snp_str);
		}		
		if(!continue_on)
			continue;
			
		id = parse_line_at_index(line, id_index);
		line_iter = find(ids.begin(), ids.end(), id);
		if(line_iter == ids.end())
			continue;			
			
		row_idx = distance(ids.begin(), line_iter);

		for(int snp = 0 ; snp < snp_vector.size(); snp++){
			snp_str = snp_vector[snp];
			SNPs(row_idx, snp) = strtod(snp_str.c_str(), NULL);
		}
	}
	
	pheno_in.close();
}
static int * get_permutation_indices(int * permutation_indices, int n_rows){
	iota (permutation_indices, permutation_indices + n_rows, 0);
	random_shuffle(permutation_indices, permutation_indices + n_rows);
//	return permutation_indices;

}	
static void calculate_pvalue(double  * pvalues, Eigen::MatrixXd syP_Sigma_P, Eigen::MatrixXd sigmaP, Eigen::MatrixXd snps_matrix){
    
	Eigen::ArrayXXd Ts = pow((syP_Sigma_P*snps_matrix).array(), 2.0)/( sigmaP*pow(snps_matrix.array(), 2.0).matrix().transpose()).array();
	//double compare_value;
	double value;
	
	for(int col = 0 ; col < Ts.cols();col++){
		double pvalue = 0.0;
		double compare_value = Ts(0, col);
		
		pvalues[col] = (Ts.col(col) >= compare_value).count()/Ts.rows();
	}

}	
static void print_gwas_help(Tcl_Interp * interp){
	Solar_Eval(interp, "help gwas");
}
extern "C" int gwaCmd(ClientData clientData, Tcl_Interp *interp,
                                         int argc,const char *argv[]){
    Eigen::initParallel();
	int n_permutations = 0;
	string plink_filename;
	vector<string> args;
	for(int arg = 1; arg < argc ; arg++){
		args.push_back(string(argv[arg]));
	}
	
	for(vector<string>::iterator it = args.begin(); it != args.end() ; it++){
		transform(it->begin(), it->end(), it->begin(), ::toupper);
		if((*it == "-NP" || *it == "--NP") && (it + 1) != args.end()){
			it++;
			n_permutations = strtol(it->c_str(), NULL, 10);
			if(n_permutations <= 0){
				RESULT_LIT("Invalid number of permutations selected.");
				return TCL_ERROR;
			}
		}else if((*it == "-H" || *it == "--H" || *it == "--HELP" || *it == "-HELP" || *it == "HELP" || *it == "H"  ) && (it + 1) != args.end()){
			print_gwas_help(interp);
			return TCL_OK;
		}else if((*it == "-PLINK" || *it == "--PLINK") && (it + 1) != args.end()){
			it++;
			plink_filename = *it;
		}else{
			RESULT_LIT("Invalid argument entered");
			print_gwas_help(interp);
			return TCL_ERROR;
		}
	}
					 
	Eigen::MatrixXd trait_vector;
	Eigen::MatrixXd eigenvectors;
	Eigen::MatrixXd aux;
	int success;				
											 
											
	if(Trait::Number_Of () == 0) {
		RESULT_LIT("Please select a trait.\n");
		return TCL_ERROR;
	}
	string current_trait = string(Trait::Name(0));
	vector<string> snp_list;
    int8_t * snp_values;
	string phenotype_filename = Phenotypes::filenames();    
	int n_subjects;
	success = generate_gwas_matrices(interp, trait_vector, &snp_values,  eigenvectors, aux, 
		 current_trait, snp_list,  plink_filename,  phenotype_filename,  n_subjects);
		 
	if(success == TCL_ERROR) return success;
	fprintf(stderr, "Running GWAS...\n");	 

	trait_vector = trait_vector.array() - trait_vector.mean();
	eigenvectors = eigenvectors.transpose().eval();
	
	trait_vector = eigenvectors*trait_vector;
//	trait_vector = trait_vector.array()/(trait_vector.norm()/sqrt(trait_vector.rows() - 1));

	// R = Hxeigenvectors.transpose()*trait_vector;

	//trait_vector = Hx_eigenvectors.transpose()*trait_vector;
	/*for(int col = 0; col < snps_matrix.cols(); col++){
		double mean_value = 0.0;
		int total_values = 0;
		for(int row = 0; row < snps_matrix.rows(); row++){
			if(snps_matrix(row, col) != -1){
				mean_value += snps_matrix(row, col)/2.0;
				total_values++;
			}
		}
		
		mean_value /= total_values;
		for(int row = 0; row < snps_matrix.rows(); row++){
			if(snps_matrix(row, col) != -1){
				snps_matrix(row, col) = snps_matrix(row, col)/2.0 - mean_value; 
			}else{
				snps_matrix(row, col) = 0.0;
			}
		}		
		
	}*/
//	snps_matrix = eigenvectors.transpose()*snps_matrix;
	
//	Eigen::MatrixXd res = trait_vector;
	Eigen::MatrixXd ZTZI = (aux.transpose()*aux).inverse();
	Eigen::VectorXd Sigma(2);
	calculate_sigma_A_and_sigma_E_fphi(trait_vector, aux, ZTZI, Sigma);
	if(Sigma(0) < 0.0 )  Sigma(0) = 0.0;
	if(Sigma(1) < 0.0 )  Sigma(1) = 0.0;
	if((Sigma(0) == 0.0 && Sigma(1) == 0.0) || (Sigma(0) != Sigma(0) || Sigma(1) != Sigma(1))){
		Sigma(0) = 1.0;
		Sigma(1) = 0.0;
	}

	Eigen::MatrixXd sigma_inverse = aux*Sigma;
	sigma_inverse = 1.0/sigma_inverse.array();
	Eigen::MatrixXd syP(trait_vector.rows(), n_permutations + 1);
	Eigen::MatrixXd sigma_P(sigma_inverse.rows(), n_permutations + 1);
	syP.col(0) = trait_vector.col(0);
	sigma_P.col(0) = sigma_inverse.col(0);
	
	
	int * indices = new int [syP.rows()];
	iota (indices, indices + n_subjects, 0);
	for(int col = 1 ; col < syP.cols(); col++){
	
	//	random_shuffle(begin(syP.col(col).data()), end(syP.col(col).data()));
		random_shuffle(indices, indices + n_subjects);
		
		for(int row = 0; row < syP.rows() ; row++){
			syP(row, col) = trait_vector(indices[row], 0);
			sigma_P(row, col) = sigma_inverse(indices[row], 0);	
				
		}
	
	}
	Eigen::MatrixXd syP_Sigma_P = (sigma_P.array()*syP.array()).matrix();
   // syP_Sigma_P = syP_Sigma_P.transpose();
   // sigma_P = sigma_P.transpose();
	string output_filename = "./" + current_trait + "/" + "gwas.out";	
	
	if(access(current_trait.c_str(), F_OK)  == -1){
		mkdir(current_trait.c_str(), 0700);
	}
	ofstream file_out(output_filename);
	file_out << "snp,pvalue\n";
	//double mean_value;
	//int total_values;
	int batch_size = GWAS_BATCH_SIZE;
	int free_snps;
	int iterations;
	if(batch_size >= snp_list.size()){
		batch_size = snp_list.size();
		iterations = 1;
		free_snps = snp_list.size();
	}else{
		iterations = ceil(float(snp_list.size())/batch_size);
		free_snps = snp_list.size() % batch_size;
	}
	int n_snps;
	double * pvalues = new double[batch_size];
	auto start = std::chrono::high_resolution_clock::now();
	vector<string>::iterator snp_iter = snp_list.begin();
	syP_Sigma_P = syP_Sigma_P.transpose();
	sigma_P = sigma_P.transpose();
	#pragma omp parallel for
	for(int iteration = 0; iteration < iterations ; iteration++){
		if(iteration < iterations - 1){
			n_snps = batch_size;
		}else{
			n_snps = free_snps;
		}
		Eigen::MatrixXd snps_matrix(n_subjects, n_snps);
		
		for(int snp = 0; snp < n_snps; snp++){
			double mean_value = 0.0;
			int total_values = 0;
			double current_value;
			for(int row = 0; row < n_subjects; row++){
				if((current_value = snp_values[(iteration*batch_size + snp)*n_subjects + row]/2.0) != 3.0/2.0){
					snps_matrix(row, snp) = current_value;
					mean_value += current_value;
					total_values++;
				}else{
					snps_matrix(row, snp) = 3.0;
				}
			}	
		
			mean_value /= total_values;
			for(int row = 0; row < n_subjects; row++){
				if(snps_matrix(row, snp) != 3.0){
					snps_matrix(row, snp) -= mean_value; 
				
				}else{
					snps_matrix(row, snp) = 0.0;
				}
			}		
		}
		snps_matrix = eigenvectors*snps_matrix;
		
		calculate_pvalue(pvalues, syP_Sigma_P,  sigma_P, snps_matrix);
		for(int snp =0; snp< n_snps; snp++, snp_iter++){
			file_out << *snp_iter << "," << pvalues[snp] << "\n";
		}	
	}
	delete [] pvalues;
	auto elapsed = std::chrono::high_resolution_clock::now() - start;
	
  auto seconds = std::chrono::duration_cast<std::chrono::duration<double>>(elapsed);	
	cout << seconds.count() << endl;
		
	delete [] indices;
	delete [] snp_values;
	return TCL_OK;			 
}
