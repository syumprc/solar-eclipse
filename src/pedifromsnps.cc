#include "solar.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <plinkio.h>
using namespace std;
static vector<const char *> read_headers(string header_filename){
	ifstream headers_in(header_filename.c_str());


	vector<const char*> headers;

	string header;


	while(headers_in >> header) headers.push_back(header.c_str());


	return headers;


}

static string to_lower(string input){

  string output;

  for(size_t idx = 0 ; idx < input.length() ; idx++) output += tolower(input[idx]);

  return output;
}
extern "C" void create_relationship_matrix(pio_file_t * input_file, const char * base_filename, double * GRM, \
 int per_chromosome, size_t n_subjects , size_t n_snps);
 
extern "C" void create_relationship_matrix_with_header(pio_file_t * input_file, double * GRM, const char * basename,
		const char * header_list[],  int n_headers, int per_chromosome, size_t n_subjects , size_t n_snps); 
 
extern "C" void write_matrix_to_file(char * output_filename, pio_file_t * input_file,  double * GRM, size_t n_subjects){

	ofstream output_file(output_filename);
	pio_sample_t * row_sample;
	pio_sample_t * col_sample;
	output_file << "IDA,IDB,KIN\n";
	for(size_t row = 0; row < n_subjects ; row++){
		row_sample = fam_get_sample(&input_file->fam_file, row);
		string row_id = row_sample->iid;
		for(size_t col = row ; col < n_subjects ; col++){
			col_sample = fam_get_sample(&input_file->fam_file, col);
			string col_id = col_sample->iid;

			if(row_id != col_id){
				output_file << row_id << "," << col_id << "," << GRM[row*n_subjects + col] << "\n";
				output_file << col_id << "," << row_id << "," << GRM[row*n_subjects + col]  << "\n";

			}else{
				output_file << row_id << "," << col_id << "," << GRM[row*n_subjects + col]  << "\n";
			}


		}
	}
	output_file.close();
}


static void print_help(Tcl_Interp * interp){
  Solar_Eval(interp, "help pedifromsnps");

}


extern "C" int pedfromsnpsCmd(ClientData clientData, Tcl_Interp *interp,
                                         int argc,const char *argv[]){

	  string args[argc];

	  for(int arg = 0 ; arg < argc ; arg++){
	      args[arg] = argv[arg];
	  }

	  string input_basename;
	  string output_basename;

	  string header_filename;
	  bool is_header = false;
	  int per_chromosome = false;
	  //Processing commandline arguments
	  for(int i = 1 ; i < argc ;i++){
	      string lower_version = to_lower(args[i]);
	      if(lower_version.compare("--help") == 0 || lower_version.compare("-help")  == 0 ||
		  lower_version.compare("-h")  == 0 || lower_version.compare("help")  == 0
		  || lower_version.compare("--h") == 0){
		  print_help(interp);
	          return EXIT_SUCCESS;
	      }else if (lower_version.compare("-i") == 0 ||
		  lower_version.compare("--i") == 0 ||
		  lower_version.compare("-input") == 0 ||
		  lower_version.compare("--input") == 0){

		  if((i+1) < argc){
		      i += 1;
		      input_basename = args[i];


		  }else{
		      printf("No input file specified.\n");
		      print_help(interp);
	              return TCL_ERROR;
		  }

	      }else if (lower_version.compare("-o") == 0 ||
		  lower_version.compare("--o") == 0 ||
		  lower_version.compare("-output") == 0 ||
		  lower_version.compare("--output") == 0){

		  if((i+1) < argc){
		      i += 1;
		      output_basename = args[i];

		  }else{

		      printf("No output file specified.\n");
		      print_help(interp);
	              return TCL_ERROR;
		  }


		  }else if ((lower_version.compare("-header") == 0 ||
				  lower_version.compare("--header") == 0 ||
				  lower_version.compare("-head") == 0 ||
				  lower_version.compare("--head") == 0) && (i + 1) < argc){


			  header_filename = args[ i + 1];
			  i++;
			  is_header = true;


	      }else if (lower_version.compare("--perchromo") == 0 || lower_version.compare("-perchromo") == 0){
			  
			  per_chromosome = 1;
			  
			  
		  }else{
			printf("Argument %s is not a valid argument.\n", args[i].c_str());
			print_help(interp);
	          return TCL_ERROR;
	      }

	  }
	  struct pio_file_t input;






	 // pio_transpose(argv[1], argv[2]);
	  pio_status_t status = pio_open(&input, input_basename.c_str());

	  if(status != PIO_OK){
	      pio_close(&input);
	      if(status == P_FAM_IO_ERROR){
		  printf("Error in loading .fam file.\n");
	          return TCL_ERROR;
	      }else if (status == P_BIM_IO_ERROR){
		  printf("Error in loading .bim file.\n");
	          return TCL_ERROR;
	      }else if (status == P_BED_IO_ERROR){
		  printf("Error in loading .bed file.\n");
	          return TCL_ERROR;
	      }
	  }

	  int version = input.bed_file.header.version;
	  if (input.bed_file.header.snp_order == BED_UNKNOWN_ORDER){
	      pio_close(&input);
	      printf("Error in the .bed snp order.\nRetry after using a different version"
		  "of plink.\n");

	      return EXIT_FAILURE;
	  }else if (input.bed_file.header.snp_order == BED_ONE_SAMPLE_PER_ROW){
	      pio_close(&input);
	      printf("In order to read efficiently the transpose of %s must be written.\n", input_basename.c_str());
	      printf("The transpose will be written to %s_trans\n", input_basename.c_str());
	      string trans_input_basename = input_basename + "_trans";
	      pio_status_t status = pio_transpose(input_basename.c_str(), trans_input_basename.c_str());
	      if(status != PIO_OK){
		  printf("Error in creating transpose.\n");
	          return TCL_ERROR;
	      }

	      status = pio_open(&input, trans_input_basename.c_str());

	      if(status != PIO_OK){
		  printf("Error in opening transpose.\n");
	          return TCL_ERROR;
	      }
	  }



	  size_t n_subjects = pio_num_samples(&input);
	  size_t n_snps = pio_num_loci(&input);


	  double * GRM = new  double[n_subjects*n_subjects];
	  memset(GRM, 0, sizeof(double)*n_subjects*n_subjects);

	  if(is_header){
		  vector<const char *>header_list = read_headers(header_filename);
		  create_relationship_matrix_with_header(&input, GRM, output_basename.c_str(), 
			header_list.data(),header_list.size(), per_chromosome,n_subjects , n_snps);
	  }else{
		  create_relationship_matrix(&input,  output_basename.c_str(), GRM,per_chromosome, n_subjects  ,  n_snps);
	  }

	 

	  pio_close(&input);
	  delete [] GRM;

	  return TCL_OK;



}

