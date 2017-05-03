#include "solar.h"
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <chrono>
#include "plinkio.h"
#include <list>
#include <unordered_map>
#include "hdf5.h"
using namespace std;
//int generate_gwas_matrices(Tcl_Interp * interp, Eigen::MatrixXd  &trait_vector, Eigen::MatrixXd  & cov_matrix,Eigen::MatrixXd & l_Hx_aux, Eigen::MatrixXd & l_aux, Eigen::MatrixXd & l_Hx_eigenvectors,\
Eigen::MatrixXd & l_eigenvectors, Eigen::MatrixXd & snps_matrix, vector<string> & snp_names, string * trait_names, int n_traits);
void calculate_eigenvectors (double * phi2, double * eigenvalues, double * eigenvectors ,int n);				
extern "C" void call_fphi(float * Y, float * Sigma, float * Z, float * hat, float * evectors,  vector<int> selected_devices,\
		int n_subjects, int n_traits, bool covariates);
extern "C" void call_gwas(float * snps, float * Y, float * Sigma, const float * evectors,  int * permutation_indices, int n_traits, \
		int n_snps, int n_permutations, int n_subjects, vector<int> selected_devices, hid_t & pvalues_out, string * trait_name);
extern "C" std::vector<int> select_devices();	
extern "C" std::vector<int> select_all_devices();	
static void print_gpu_gwas_help(Tcl_Interp * interp){
	Solar_Eval(interp, "help gpu_gwas");	
}

static int generate_gpu_gwas_matrices(Tcl_Interp * interp, float ** trait_values, float ** snp_values, float ** evectors, float ** aux, 
		vector<string> trait_names, vector<string> & snp_list, string plink_filename, string phenotype_filename, int & n_subjects){
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
	
	for(vector<string>::iterator trait = trait_names.begin(); trait != trait_names.end() ; trait++){
		file->setup(trait->c_str(), &errmsg);				
		if(errmsg){
			RESULT_LIT(errmsg);
			return TCL_ERROR;
		}
	}
	unordered_map<string , list<double> > trait_map;
	list<double> row_values(trait_names.size());
	list<string>::iterator id_iter = phenotype_ids.begin();
	list<double>::iterator row_iter;
	bool skip_row;
	list<string> bad_ids;
	int n_traits = trait_names.size();
	
	while (0 != (file_data = file->get (&errmsg))){
		row_iter = row_values.begin();
		skip_row = false;
		for(int idx = 0; idx < n_traits; idx++, row_iter++){
			if(StringCmp(file_data[idx],0, 0)){
				*row_iter = strtod((const char*)file_data[idx], NULL);
			}else{
				skip_row = true;
				break;
			}
		}
		if(skip_row){
			
			bad_ids.push_back(*id_iter++);	
			continue;
		}
		
		trait_map[*id_iter] = row_values;
		id_iter++;
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
	float * temp_Y = new float[n_subjects*n_traits];
	int trait_idx;
	int id_index = 0;
	for(id_iter = ids.begin(); id_iter != ids.end(); id_iter++, id_index++){
		row_values = trait_map[*id_iter];
		trait_idx = 0;
		for(row_iter = row_values.begin(); row_iter != row_values.end(); row_iter++, trait_idx++){
			temp_Y[trait_idx*ids.size() + id_index] = *row_iter;
		}
		
	}
	trait_map.clear();
	*trait_values = temp_Y;
	
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
	
	float * temp_evectors = new float[n_subjects*n_subjects];
	float * temp_Z = new float[n_subjects*2];
	
	for(int row = 0; row < n_subjects; row++){
		temp_Z[row] = 1.f;
		temp_Z[n_subjects + row] = (float)eigenvalues[row];
		for(int col = 0; col < n_subjects; col++){
			temp_evectors[col*n_subjects + row] = (float)d_eigenvectors[col*n_subjects + row];
		}
	}
	delete [] eigenvalues;
	delete [] d_eigenvectors;
	*evectors = temp_evectors;
	*aux = temp_Z;
	n_snps = snp_list.size();	
	if(plink_filename != ""){
		snp_t * snp_buffer = new snp_t[sample_size];
		float * temp_snps = new float[n_snps*ids.size()];
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
		float * temp_snps = new float[n_snps*ids.size()];
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

static vector<string>  get_headers(string header_filename){
	vector<string> header_list;
	string header;
	ifstream header_in(header_filename.c_str());

	while(header_in >> header)
		header_list.push_back(header);

	header_in.close();

	return header_list;

}
/*
static float * convert_to_pointer(Eigen::MatrixXd matrix){
	float * pointer = new float[matrix.cols()*matrix.rows()];
	for(int col = 0 ; col < matrix.cols() ;col++){
		for(int row = 0; row < matrix.rows(); row++){
			pointer[col*matrix.rows() + row] = matrix(row, col);
		}
	}
	return pointer;
}*/
static int * create_permutation_indices(int n_subjects, int n_permutations){
	int * permutation_indices = new int[n_subjects*n_permutations];
	for(int permutation = 0; permutation < n_permutations ; permutation++){
		iota (permutation_indices + permutation*n_subjects, permutation_indices + permutation*n_subjects + n_subjects, 0);
		random_shuffle (permutation_indices+ permutation*n_subjects, permutation_indices + permutation*n_subjects + n_subjects);	
	}		
	return permutation_indices;	
}
static vector<string> get_snp_names(){
	ifstream snps_in("snp.geno-list");
	string snp_name;
	vector<string> snp_list;
	while(snps_in >> snp_name) snp_list.push_back(snp_name);
	
	snps_in.close();
	
	return snp_list;
}

 
static void initialize_output_file(const char * filename, hid_t & pvalues_out, vector<string> snp_names, int n_snps){
	
	pvalues_out  = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	hsize_t dim[1] = {n_snps};
	
	hid_t datatype = H5Tcopy (H5T_C_S1);
	H5Tset_size(datatype, H5T_VARIABLE);
	hid_t dataspace = H5Screate_simple (1, dim, NULL);

	hid_t plist_id  = H5Pcreate (H5P_DATASET_CREATE);

    H5Pset_chunk (plist_id, 1, dim);
   
    H5Pset_deflate (plist_id, 9);

    hid_t dataset_id = H5Dcreate2 (pvalues_out, "/snps", datatype,
                            dataspace, H5P_DEFAULT, plist_id, H5P_DEFAULT);


    H5Dwrite (dataset_id, datatype, dataspace, H5S_ALL, H5P_DEFAULT, snp_names.data());
   

    H5Sclose (dataspace);
    H5Dclose (dataset_id);
    H5Pclose (plist_id);
	
}		
extern "C" int gpuGWASCmd(ClientData clientData, Tcl_Interp *interp, 
						int argc , const char * argv[]){
	
	int n_permutations;
	string header_name;
	string output_filename;
	string plink_filename;
	vector<string> args;
	bool select_all = false;
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
		}else if((*it == "-HEAD" || *it == "--HEADER" || *it == "--HEAD" || *it == "-HEADER"  ) && (it + 1) != args.end()){
			it++;
			header_name = *it;
		}else if((*it == "-H" || *it == "--H" || *it == "--HELP" || *it == "-HELP" || *it == "HELP" || *it == "H"  ) && (it + 1) != args.end()){
			print_gpu_gwas_help(interp);
		}else if((*it == "-O" || *it == "--O" || *it == "--OUT" || *it == "-OUT" || *it == "--OUTPUT" || *it == "-OUTPUT" ) && (it + 1) != args.end()){
			it++;
			output_filename = *it;
		}else if((*it == "-PLINK" || *it == "--PLINK") && (it + 1) != args.end()){
			it++;
			plink_filename = *it;
		}else if (*it == "-ALL" || *it == "--ALL" || *it == "--A" || *it == "-A") {
			select_all = true;
			
		}else{
			RESULT_LIT("Invalid argument entered");
			return TCL_ERROR;
		}
	}
	vector<string> headers = get_headers(header_name);
	vector<string> snp_list;
	string phenotype_filename =  Phenotypes::filenames();
	float * Y, * evectors, * Z, * hat, * snps; 	
	int n_subjects, n_traits, n_snps;		
	n_traits = headers.size();

	int success = generate_gpu_gwas_matrices(interp, &Y, &snps, &evectors, &Z, 
		headers, snp_list, plink_filename,phenotype_filename, n_subjects);			
	if(success == TCL_ERROR) return success;		
//	int success = generate_gpu_gwas_matrices(interp, &Y, &snps, &evectors, &Z, 
//		headers, snp_list,  plink_filename,  phenotype_filename, n_subjects);
 //   if(success == TCL_ERROR) return success;
    n_snps = snp_list.size();
    if(n_subjects == 0){
		RESULT_LIT("No subjects were found.");
		return TCL_ERROR;	
	} 
    if(n_snps == 0){
		RESULT_LIT("No snps were found.");
		return TCL_ERROR;	
	} 	
    if(n_traits == 0){
		RESULT_LIT("No traits were found.");
		return TCL_ERROR;	
	}
	vector<int> selected_devices;
	if(!select_all){	
		selected_devices = select_devices();
	}else{
		selected_devices = select_all_devices();
	}
 
		
 
    
    bool covariates = false;
       

    float * Sigma = new float[n_traits*n_subjects]; 
	
	auto start = std::chrono::high_resolution_clock::now();
	call_fphi(Y, Sigma, Z,  hat,  evectors, selected_devices,\
		 n_subjects,  n_traits,  covariates);
	if(covariates) delete [] hat;
	delete [] Z;
	hid_t pvalues_out;
	initialize_output_file(output_filename.c_str(), pvalues_out, snp_list, n_snps);	    
    int * permutation_indices = create_permutation_indices(n_subjects ,n_permutations);	
    	
	call_gwas(snps, Y, Sigma, evectors, permutation_indices, n_traits, \
		 n_snps,  n_permutations,  n_subjects, selected_devices, pvalues_out, headers.data());
	H5Fclose(pvalues_out);
	 auto elapsed = std::chrono::high_resolution_clock::now() - start;
	
  auto seconds = std::chrono::duration_cast<std::chrono::duration<double>>(elapsed);	
	cout << seconds.count() << endl;		 
		 
	delete [] Y;
	delete [] Sigma;
	delete [] evectors;
	delete [] permutation_indices;
	//delete [] snps;

	return TCL_OK;	
}		
