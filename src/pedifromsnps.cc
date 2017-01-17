#include<plinkio.h>
#include<string>
#include<fstream>
#include<stdio.h>
#include<iostream>
#include <algorithm>
#include <vector>
#include <iterator>
#include "solar.h"
using namespace std;
#define SQUARED(value) ((value)*(value))


static string to_lower(string input){

  string output;

  for(size_t idx = 0 ; idx < input.length() ; idx++) output += tolower(input[idx]);

  return output;
}



inline static long unsigned int decrease (long unsigned int i) { return --i; }

static void create_relationship_matrix(pio_file_t * input_file, vector <vector <long double> > & GRM, size_t n_subjects , size_t n_snps){

	snp_t * buffer = new snp_t[n_subjects];
	// * int_buffer = new unsigned int[n_subjects];




	 long unsigned  int * n_snps_matrix= new long unsigned int[n_subjects*n_subjects];

	fill(n_snps_matrix, n_snps_matrix + n_subjects*n_subjects, n_snps);


	pio_locus_t * locus;
	pio_status_t status;


	double freq;


	 unsigned int homo_count;
	 unsigned int heter_count;
	 unsigned int subject_count;

	size_t row, col, current_subject;




	vector<unsigned char> vector_buffer(n_subjects);
	for(size_t snp = 0 ; snp < n_snps; ++snp){



		status = pio_next_row(input_file,buffer);
		if(status != PIO_OK){
			for(size_t row = 0 ; row < n_subjects ; ++row) \
					for(size_t col = row;  col < n_subjects; ++col) n_snps_matrix[row*n_subjects + col]--;
			continue;
		}

		locus = bim_get_locus(&input_file->bim_file, snp);
		if(strcmp(locus->allele1, "0") == 0 || strcmp(locus->allele2, "0")== 0){
			for(size_t row = 0 ; row < n_subjects ; ++row) \
					for(size_t col = row;  col < n_subjects; ++col) n_snps_matrix[row*n_subjects + col]--;
			continue;
		}
		freq = 0.0;


		subject_count = n_subjects;
		homo_count = 0;
		heter_count = 0;


		copy(buffer, buffer + n_subjects, vector_buffer.begin());





		auto it_buffer = vector_buffer.begin();

		for(size_t subject = 0 ; subject < n_subjects ; ++subject){

			switch(*it_buffer++){
			case 1:
				heter_count++;
				break;
			case 2:
				homo_count++;
				break;
			case 3:
				subject_count--;

				for(size_t col = subject; col < n_subjects ; ++col) n_snps_matrix[subject*n_subjects + col]--;

				for(size_t row = 0; row < subject; ++row) n_snps_matrix[row*n_subjects + subject]--;




				break;

			}




		}

		freq = (double)(homo_count*2.0 + heter_count)/(double)(2.0*subject_count);



		auto row_value = vector_buffer.begin();
		for(size_t row = 0; row < n_subjects; ++row, ++row_value){
			if(*row_value == 3) continue;
			GRM[row][row] += SQUARED((*row_value - 2.0*freq))/(freq*2.0*(1.0 - freq));


			auto col_value = vector_buffer.end() - 1;
			double factor = *row_value - 2.0*freq;
			for(size_t col = n_subjects - 1 ; col > row ; --col, --col_value) GRM[row][col] += factor*((*col_value != 3 ? *col_value : 2.0*freq)  - \
							2.0*freq)/(freq*(1.0 - freq)*2.0);





		}



	}
	delete [] buffer;


	for(size_t row = 0 ; row < n_subjects ; ++row)\
			for(size_t col = row;  col < n_subjects; ++col) GRM[row][col] /= (long double)n_snps_matrix[row*n_subjects + col];




	delete [] n_snps_matrix;

}

static void create_relationship_matrix_with_header(pio_file_t * input_file, vector <vector <long double> > & GRM,
		vector<string> header_list, size_t n_subjects , size_t n_snps){

	snp_t * buffer = new snp_t[n_subjects];
	// * int_buffer = new unsigned int[n_subjects];




	 long unsigned  int * n_snps_matrix= new long unsigned int[n_subjects*n_subjects];

	fill(n_snps_matrix, n_snps_matrix + n_subjects*n_subjects, n_snps);


	pio_locus_t * locus;
	pio_status_t status;


	double freq;

	 unsigned int homo_count;
	 unsigned int heter_count;
	 unsigned int subject_count;
	size_t row, col, current_subject;




	vector<unsigned char> vector_buffer(n_subjects);
	for(size_t snp = 0 ; snp < n_snps; ++snp){



		status = pio_next_row(input_file,buffer);
		if(status != PIO_OK){
			for(size_t row = 0 ; row < n_subjects ; ++row) \
					for(size_t col = row;  col < n_subjects; ++col) n_snps_matrix[row*n_subjects + col]--;
			continue;
		}

		locus = bim_get_locus(&input_file->bim_file, snp);
		if(strcmp(locus->allele1, "0") == 0 || strcmp(locus->allele2, "0")== 0){
			transform (n_snps_matrix, n_snps_matrix + n_subjects*n_subjects, n_snps_matrix, decrease);
			continue;
		}
		freq = 0.0;

		string header_name(locus->name);

		if(find(header_list.begin(), header_list.end(), header_name) == header_list.end()){

			for(size_t row = 0 ; row < n_subjects ; ++row) \
					for(size_t col = row;  col < n_subjects; ++col) n_snps_matrix[row*n_subjects + col]--;
			continue;
		}


		subject_count = n_subjects;
		homo_count = 0;
		heter_count = 0;


		copy(buffer, buffer + n_subjects, vector_buffer.begin());





		auto it_buffer = vector_buffer.begin();

		for(size_t subject = 0 ; subject < n_subjects ; ++subject){

			switch(*it_buffer++){
			case 1:
				heter_count++;
				break;
			case 2:
				homo_count++;
				break;
			case 3:
				subject_count--;

				for(size_t col = subject; col < n_subjects ; ++col) n_snps_matrix[subject*n_subjects + col]--;

				for(size_t row = 0; row < subject; ++row) n_snps_matrix[row*n_subjects + subject]--;




				break;

			}




		}

		freq = (double)(homo_count*2.0 + heter_count)/(double)(2.0*subject_count);



		auto row_value = vector_buffer.begin();
		for(size_t row = 0; row < n_subjects; ++row, ++row_value){
			if(*row_value == 3) continue;
			GRM[row][row] += SQUARED((*row_value - 2.0*freq))/(freq*2.0*(1.0 - freq));


			auto col_value = vector_buffer.end() - 1;
			double factor = *row_value - 2.0*freq;
			for(size_t col = n_subjects - 1 ; col > row ; --col, --col_value) GRM[row][col] += factor*((*col_value != 3 ? *col_value : 2.0*freq)  - \
							2.0*freq)/(freq*(1.0 - freq)*2.0);





		}



	}
	delete [] buffer;


	for(size_t row = 0 ; row < n_subjects ; ++row) \
			for(size_t col = row;  col < n_subjects; ++col) GRM[row][col] /= (long double)n_snps_matrix[row*n_subjects + col];




	delete [] n_snps_matrix;

}

static vector<string> read_headers(string header_filename){
	ifstream headers_in(header_filename.c_str());


	vector<string> headers;

	string header;


	while(headers_in >> header) headers.push_back(header);


	return headers;


}

static void write_matrix_to_file(string output_filename, pio_file_t * input_file, vector< vector<long double> > GRM, size_t n_subjects){

	ofstream output_file(output_filename.c_str());
	pio_sample_t * row_sample;
	pio_sample_t * col_sample;
	output_file << "IDA,IDB,KIN\n";
	for(size_t row = 0; row < n_subjects ; row++){
		row_sample = fam_get_sample(&input_file->fam_file, row);
		string row_id = row_sample->iid;
		for(size_t col = row ; col < n_subjects ; col++){
			col_sample = fam_get_sample(&input_file->fam_file, col);
			string col_id = col_sample->iid;

			if(row_id.compare(col_id) == 1){
				output_file << row_id << "," << col_id << "," << (long double)GRM[row][col] << "\n";
				output_file << col_id << "," << row_id << "," << (long double)GRM[row][col] << "\n";

			}else{
				output_file << row_id << "," << col_id << "," <<  (long double)GRM[row][col] << "\n";
			}


		}
	}
	output_file.close();
}
static void print_help(){
  printf("Use: help pedifromsnps\n");
}
extern "C" int Runpedfromsnps(ClientData clientData, Tcl_Interp *interp,
                                         int argc,const char *argv[]){
	  if(argc != 5 && argc != 6 && argc != 8 && argc != 2){
	      printf("Invalid number of arguments\n");
	      print_help();
	      return TCL_OK;
	  }


	  string args[argc];

	  for(int arg = 0 ; arg < argc ; arg++){
	      args[arg] = argv[arg];
	  }

	  string input_basename;
	  string output_basename;

	  string header_filename;
	  bool is_header = false;
	  //Processing commandline arguments
	  for(int i = 1 ; i < argc ;i++){
	      string lower_version = to_lower(args[i]);
	      if(lower_version.compare("--help") == 0 || lower_version.compare("-help")  == 0 ||
		  lower_version.compare("-h")  == 0 || lower_version.compare("help")  == 0
		  || lower_version.compare("--h") == 0){
		  print_help();
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
		      print_help();
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
		      print_help();
	              return TCL_ERROR;
		  }


		  }else if ((lower_version.compare("-header") == 0 ||
				  lower_version.compare("--header") == 0 ||
				  lower_version.compare("-head") == 0 ||
				  lower_version.compare("--head") == 0) && (i + 1) < argc){


			  header_filename = args[ i + 1];
			  i++;
			  is_header = true;


	      }else{
		  printf("Argument %s is not a valid argument.\n", args[i].c_str());
		  print_help();
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


	  vector< vector<long double> > GRM;

	  for(size_t i = 0 ; i < n_subjects; ++i){
		  vector <long double> GRM_row(n_subjects);
		  GRM.push_back(GRM_row);

	  }

	  if(is_header){
		  vector<string>header_list = read_headers(header_filename);
		  create_relationship_matrix_with_header(&input, GRM,
			header_list, n_subjects , n_snps);
	  }else{
		  create_relationship_matrix(&input, GRM,  n_subjects ,  n_snps);
	  }

	  write_matrix_to_file(output_basename, &input, GRM, n_subjects);

	  pio_close(&input);
	 // delete [] GRM;

	  return TCL_OK;



}


#undef SQUARED

