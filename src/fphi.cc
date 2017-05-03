#include <Eigen/Dense>
#include <vector>
#include <string>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include "solar.h"
#include <chrono>

using namespace std;


int create_eigenvectors_and_auxiliary_matrix(Tcl_Interp * interp, Eigen::MatrixXd & aux, Eigen::MatrixXd  &eigenvectors);
vector<string> fetch_ids();	
int create_trait_and_covariate_matrix(Tcl_Interp * interp, Eigen::MatrixXd & trait_matrix, Eigen::MatrixXd & covariate_matrix, \
									vector<string> unique_cov_list, vector<string> cov_list, int n_traits);
									
 int  create_id_list(string phenotype_filename,  vector<string> trts,  vector<string> covs, int n_trts, int n_covs);
extern "C" void symeig_ (int*, double*, double*, double*, double*, int*);
extern "C" void cdfchi_ (int*, double*, double*, double*, double*,
			 int*, double*);
			 
extern "C" double tdist_(double * X, double * P);			 
			 
static double chicdf(double chi, double df){
	double p, q, bound;
	int status = 0;
	int which = 1;
	
	
	cdfchi_ (&which, &p, &q, &chi, &df, &status, &bound);
	
	return q/2.0;
}	
static void calculate_eigenvectors (double * phi2, double * eigenvalues, double * eigenvectors ,int n)
{
    double* e =  new double[n];
    memset(e, 0, sizeof(double)*n);
    int * info = new int;
    *info  = 0;
    symeig_(&n, phi2, eigenvalues, e, eigenvectors, info);
    delete [] e;
    delete [] info;
}





extern "C" void ibdid2phi_ (int* ibdid1, int* ibdid2, double *phi2);
extern "C" void cdfnor_ (int* which, double* p, double* q, double* x, double* mean, double* sd, int* status, double* bound);

double Zscore_to_pvalue(double X){
	
	double mean = 0;
	int which = 1;
	double sd = 1.0;
	double q;
	double p;
	int status;
	double bound;
	
	cdfnor_(&which, &p, &q, &X, &mean, &sd, &status, &bound);
	
	if(p < q)
		return p;
	else
		return q;
		
}
	



static inline bool return_isspace (char  x){
	if(std::isspace(x))
		return true;
    else
	   return false;

}

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
	
void calculate_sigma_A_and_sigma_E_fphi(Eigen::MatrixXd Y, Eigen::MatrixXd Z, Eigen::MatrixXd ZTZI, Eigen::VectorXd & Sigma){

	Eigen::VectorXd F = pow(Y.array(), 2).matrix();
	double F_sum = F.array().mean();
	
	double score =  1.0/F_sum * (Z.col(1).array()*((F.array()/F_sum) - 1.0)).sum();
	if (score != score) {
		Sigma(0) = 0.0;
		Sigma(1) = 0.0;
		return;
	}
	if(score <= 0.0){
		Sigma(0) = 0.0;
		Sigma(1) = 0.0;
		return;
	 }
	Eigen::VectorXd theta = ZTZI*Z.transpose()*F;

	if(theta(0) < 0.f) theta(0) = 0.f;
	if(theta(1) < 0.f) theta(1) = 0.f;
	Eigen::VectorXd weights = Z*theta;

	Eigen::MatrixXd weight_matrix  = pow(weights.array(), -2).matrix().asDiagonal();

	Eigen::MatrixXd ZWZ = Z.transpose()*weight_matrix*Z;
	
	Sigma = ZWZ.fullPivLu().solve(Z.transpose()*weight_matrix*F);

}
	
void calculate_h2r(Eigen::MatrixXd Y, Eigen::MatrixXd Z, Eigen::MatrixXd ZTZI, float & h2r , float & score, float & SE){

	Eigen::VectorXd F = pow(Y.array(), 2).matrix();
	double F_sum = F.array().mean();
	
	score =  1.0/F_sum * (Z.col(1).array()*((F.array()/F_sum) - 1.0)).sum();

	Eigen::VectorXd theta = ZTZI*Z.transpose()*F;

	if(theta(0) < 0.f) theta(0) = 0.f;
	if(theta(1) < 0.f) theta(1) = 0.f;
	
	if (isnan(score)) {
		score = 0.f;
		h2r = 0.f;
		SE = 0.f;
		return;
	}
	if(score <= 0.0){
		h2r = 0.0;
		SE = 0.f;
		 return;
	 }	
	Eigen::VectorXd weights = Z*theta;


	Eigen::MatrixXd weight_matrix  = pow(weights.array(), -2).matrix().asDiagonal();

	Eigen::MatrixXd ZWZ = Z.transpose()*weight_matrix*Z;
	
	Eigen::VectorXd Sigma = ZWZ.fullPivLu().solve(Z.transpose()*weight_matrix*F);

	if(Sigma(0) < 0.f) Sigma(0) = 0.f;
	if(Sigma(1) < 0.f) Sigma(1) = 0.f;
	h2r = Sigma(1)/(Sigma(1) + Sigma(0));
	weights = Z*Sigma;
	
	if(isnan(h2r)) {
		h2r = 0.f;
		SE = 0.f;
	}


	 weights = pow((Z*Sigma).array(), -2).matrix();
	 double a = weights.sum();
	 
	 double b = (Z.col(1).array()*weights.array()).sum();
	 
	 double c = (pow(Z.col(1).array(), 2)*weights.array()).sum();
	 
	 double det = a*c - pow(b, 2.0);
	 
	 double G = Sigma(1)/pow(Sigma(1) + Sigma(0), 2);
	 double E = Sigma(0)/pow(Sigma(1) + Sigma(0), 2);
	 
	 double var =2.0*(pow(G, 2.0)*c + 2.0*(G*E)*b + pow(E, 2.0)*a)/det;
	 SE = sqrt(var);
	 


}
static void print_fphi_help(Tcl_Interp * interp){
	Solar_Eval(interp, "help fphi");
}

int generate_fphi_matrices(Tcl_Interp * interp, Eigen::MatrixXd  & l_trait_vector, Eigen::MatrixXd  & l_cov_matrix, Eigen::MatrixXd & l_aux, Eigen::MatrixXd & l_eigenvectors,\
				vector<string> trait_names, int n_traits);				
						
extern "C" int runfphiCmd(ClientData clientData, Tcl_Interp * interp,
								int argc, const char * argv[]){
									
	Eigen::MatrixXd trait_vector;
	Eigen::MatrixXd cov_matrix;
	Eigen::MatrixXd eigenvectors;
	Eigen::MatrixXd aux;
	

	int success;				
	if(argc > 1){
		print_fphi_help(interp);
		return TCL_OK;
	}
	
	if(Trait::Number_Of() == 0){
		RESULT_LIT("No trait has been selected");
		return TCL_ERROR;
	}
	
	vector<string> trait_list;
	trait_list.push_back(string(Trait::Name(0)));
	success = generate_fphi_matrices(interp, trait_vector, cov_matrix, aux, eigenvectors,trait_list , 1);
	
	if(success == TCL_ERROR) return success;
	trait_vector = trait_vector.array() - trait_vector.mean();
	trait_vector = trait_vector.array()/(trait_vector.norm()/sqrt(trait_vector.rows() - 1));
	trait_vector = eigenvectors.transpose()*trait_vector;
	
	Eigen::VectorXd R = trait_vector;
	if(cov_matrix.cols() != 0){
	
		for(int col = 0; col < cov_matrix.cols() - 1; col++){
			cov_matrix.col(col) =eigenvectors.transpose()*cov_matrix.col(col);
		}		

		//cov_matrix  = eigenvectors.transpose()*cov_matrix; 
		Eigen::MatrixXd beta = (cov_matrix.transpose()*cov_matrix).inverse()*cov_matrix.transpose()*trait_vector;
		R  = trait_vector - cov_matrix*beta;
	}
//	R = eigenvectors.transpose()*R;
	Eigen::MatrixXd ZTZI = (aux.transpose()*aux).inverse();
	float h2r;
	float score;
	float SE;
	calculate_h2r(R, aux,ZTZI, h2r , score, SE);

	
	string results = to_string(score) + " " + to_string(h2r) + " " + to_string(SE);
	
	RESULT_BUF(results.c_str());
	

	return TCL_OK;
	
}	
/*
	
extern "C" int Runcfphicmd(ClientData clientData, Tcl_Interp *interp,
                                         int argc,const char *argv[])
{

	vector<string> args;
 // std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

	for(int arg = 1; arg < argc; arg++) args.push_back(string(argv[arg]));
	Eigen::MatrixXd Z;
	Eigen::MatrixXd Y;
	Eigen::MatrixXd X;
	for(vector<string>::iterator it = args.begin() ;  it != args.end(); it++){
		transform((*it).begin(), (*it).end(), (*it).begin(), ::toupper);
		if((*it).compare("-H") == 0 || (*it).compare("--H") == 0 ||
				(*it).compare("-HELP") == 0 || (*it).compare("--HELP") == 0 ){
		
			print_cfphi_help(interp);
			return TCL_OK;
		}else{
			fprintf(stderr, "An invalid argument was entered.\n");
			print_cfphi_help(interp);
			return TCL_ERROR;
		}
	}
	
	if(Trait::Number_Of () == 0) {
		printf("Please select a trait.\n");
		return TCL_ERROR;
	}
	int success;
	string current_trait = string(Trait::Name(0));
	
	transform(current_trait.begin(), current_trait.end(), current_trait.begin(), ::toupper); 
	if(current_trait != fphi_trait || eigenvectors_matrix_index == -1 || trait_index == -1) success = create_trait_vector(interp) ;
	
	if(success == TCL_ERROR) return success;
	
	if(reload_covs || cov_matrix_index == -1) success = create_cov_matrix(interp);
	
	if(success == TCL_ERROR) return success;
	
	Z = *MathMatrixGet (CONVERT_INDEX(auxiliary_matrix_index), "FPHI", interp);
	
	Y = *MathMatrixGet (CONVERT_INDEX(trait_index), "FPHI", interp);
	if(cov_matrix_index != -1) X = *MathMatrixGet (CONVERT_INDEX(cov_matrix_index), "FPHI", interp);
		
	Eigen::MatrixXd R;
	if (X.cols() != 0) {
		Eigen::MatrixXd beta = (X.transpose()*X).ldlt().solve(X.transpose()*Y);
	
		R  = Y - X*beta;
	}else{
		R = Y;
	}
	
	
	if(R.norm()/sqrt(R.rows() - 1.0) < 0.000001) R = Eigen::ArrayXd::Zero(R.rows()).matrix();
	
		
	float h2r = 0.f;
	float score = 0.f;
	Eigen::MatrixXd ZTZI = (Z.transpose()*Z).inverse();
	calculate_h2r(R, Z, ZTZI, h2r, score);
	stringstream ss;
	ss  << score << " " << h2r;
	
	printf("Score : %f\n", score);
	printf("H2r : %f\n", h2r);
	
	RESULT_BUF(ss.str().c_str());
	
//	auto elapsed = std::chrono::high_resolution_clock::now() - start;
//	std::chrono::duration<double> seconds = std::chrono::duration_cast<std::chrono::duration<double>>(elapsed);
	

	return TCL_OK;
}
*/

