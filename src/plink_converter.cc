#include "solar.h"
#include "safelib.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <bitset>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "plinkio/status.h"
#include <cmath>
#include "plinkio/plinkio.h"


using namespace std;

struct alleles{
  char  * allele1;
  char  * allele2;
};




static void print_help(){
  printf("Use: help plink_converter\n");
}

static string to_lower(string input){

  string output;

  for(size_t idx = 0 ; idx < input.length() ; idx++){
      output += tolower(input[idx]);
  }
  return output;
}




static  string convert_to_text(int val){
  ostringstream ss;
  ss << val;
  return ss.str();
}

static int read_and_write_to_csv( pio_file_t  * input, string output_basename,
			  size_t col_size, size_t max_per_file, bool bin_format){

  size_t row_size = input->bed_file.header.num_samples;
  int num_files;
  if (max_per_file == 0){
      max_per_file = col_size;
      num_files = 1;
  }else{
      num_files = std::ceil(float(col_size)/float(max_per_file));
  }
  alleles current_allele;


  snp_t buffer[row_size];
  pio_locus_t * locus;
  ofstream  oss[num_files];


  if (num_files == 1){


      string output_filename = output_basename + ".csv";
      oss[0].open(output_filename.c_str(), ios::out);


  }else{

      for(unsigned int i = 0 ; i < num_files; i++){
	  string output_filename = output_basename + "_" + convert_to_text(i) + ".csv";
	  oss[i].open(output_filename.c_str(), ios::out);
      }


  }




  oss[0] << "id,fid,sex,phenotype,";

  int file = 0;
  bitset<8> bits;

  for(unsigned int i = 0 ; i < col_size ; i++){

      if((max_per_file*(file + 1)) == i){
	  file++;

	  oss[file] << "id,fid,sex,phenotype,";
      }

      locus = bim_get_locus(&input->bim_file, i);


      if(i == ((file + 1)*(max_per_file) - 1) || (i == col_size - 1)){
	  oss[file] << "snp_" << locus->name << "\n";
      }else{
	  oss[file] << "snp_" << locus->name << ",";
      }

  }

  pio_sample_t * sample;
  for(size_t row = 0 ; row < row_size ; row++){


      sample = fam_get_sample(&input->fam_file, row);
      file = 0;
      oss[file] << sample->iid << "," << sample->fid << "," << sample->sex << "," << sample->phenotype << ",";
      for(size_t col = 0 ; col < col_size ; col++){
	  pio_status_t status = pio_next_row(input, buffer);


	   if((max_per_file*(file + 1)) == col){
	       file++;
	       sample = fam_get_sample(&input->fam_file, row);
	       oss[file] << sample->iid << "," << sample->fid << "," << sample->sex << "," << sample->phenotype << ",";
	   }
	   bits = buffer[row];
	   unsigned long value = bits.to_ulong();
	   bits.reset();
	   locus = bim_get_locus(&input->bim_file, col);
    	   alleles current;
    	   current.allele1 = locus->allele1;
    	   current.allele2 = locus->allele2;
    	   if((string(current.allele1).compare("0") == 0
    	       || string(current.allele2).compare("0") == 0)  && value == 1){
    	       value = 3;
    	   }else if((string(current.allele1).compare("0") == 0 )
    		  && value == 0){
    		  value = 3;
    	   }else if((string(current.allele2).compare("0") == 0 )
 	     && value == 2){
    	       value = 3;
 	   }


    	   if(value == 0){

    	       if(bin_format == true){
    		   oss[file] << value;
    	       }else{
    		   oss[file] << current.allele1 << current.allele1;
    	       }
    	   }else if (value == 1){
    	       if(bin_format == true){
    		   oss[file] << value;
    	       }else{
    		   oss[file] << current.allele1 << current.allele2;
    	       }
    	   }else if (value == 2){
    	       if(bin_format == true){
    		   oss[file] << value;
    	       }else{
    		   oss[file] << current.allele2 << current.allele2;
    	       }

    	   }else{
    	       oss[file] << "";
    	   }
    	   if(col == ((file + 1)*(max_per_file) - 1) || (col == col_size - 1)){
    	       oss[file] << "\n";
    	   }else{
    	       oss[file] << ",";
    	   }
      }
      pio_reset_row(input);
  }

  for(int f = 0 ; f < num_files ; f++){
      oss[f].close();
  }

  pio_close(input);

  return TCL_OK;

}


extern "C" int RunPlinkConverter(ClientData clientData, Tcl_Interp *interp,
                                         int argc,const char *argv[])
{

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
  bool bin_format = false;
  size_t max_per_file = 0;

  //Processing commandline arguments
  for(int i = 1 ; i < argc ;i++){
      string lower_version = to_lower(args[i]);
      if(lower_version.compare("--help") == 0 || lower_version.compare("-help")  == 0 ||
	  lower_version.compare("-h")  == 0 || lower_version.compare("help")  == 0
	  || lower_version.compare("--h") == 0){
	  print_help();
          return TCL_OK;
      }else if(lower_version.compare("--bin") == 0 ||
	  lower_version.compare("-bin") == 0 ||
	  lower_version.compare("-b") == 0  ||
	  lower_version.compare("--b") == 0){

	  bin_format = true;

      }else if(lower_version.compare("--max") == 0 ||
	  lower_version.compare("-max") == 0 ||
	  lower_version.compare("--m") == 0||
	  lower_version.compare("-m") == 0 ){
	  if((i+1) < argc){
	      i += 1;
	      max_per_file = strtol(args[i].c_str(), NULL, 10);

	      if(max_per_file == 0){
		  printf("Invalid entry for max option or zero was entered.\n");
		  print_help();
                  return TCL_ERROR;
	      }

	  }else{
	      printf("Invalid entry for max option.\n");
	      print_help();
              return TCL_ERROR;
	  }
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

      return TCL_ERROR;
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
  size_t col_size = input.bed_file.header.num_loci;
  return read_and_write_to_csv(&input, output_basename,
  			 col_size, max_per_file,  bin_format);




}
