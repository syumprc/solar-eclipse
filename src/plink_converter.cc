#include "solar.h"
#include "safelib.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <stdio.h>
#include "plinkio.h"
#include <algorithm>

using namespace std;


static void print_help(Tcl_Interp * interp){
  Solar_Eval(interp,"help plink_converter\n");
}


extern "C" void read_and_write_to_csv(struct pio_file_t  * input, const char* output_basename,
			  size_t col_size, size_t max_per_file, int bin_format, int per_chromo);


extern "C" int RunPlinkConverter(ClientData clientData, Tcl_Interp *interp,
                                         int argc,const char *argv[])
{




  string args[argc];

  for(int arg = 0 ; arg < argc ; arg++){
      args[arg] = argv[arg];
  }

  string input_basename;
  string output_basename;
  bool bin_format = false;
  int per_chromo = 0;
  size_t max_per_file = 0;

  //Processing commandline arguments
  for(int i = 1 ; i < argc ;i++){
      string lower_version = args[i];
      transform(lower_version.begin(), lower_version.end(), lower_version.begin(), ::tolower);
      if(lower_version.compare("--help") == 0 || lower_version.compare("-help")  == 0 ||
	  lower_version.compare("-h")  == 0 || lower_version.compare("help")  == 0
	  || lower_version.compare("--h") == 0){
	  print_help(interp);
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
		  print_help(interp);
                  return TCL_ERROR;
			  }
      }else{
	      printf("Invalid entry for max option.\n");
	      print_help(interp);
              return TCL_ERROR;
		  }
	  }else if(lower_version.compare("--per-chromo") == 0 ||
	  lower_version.compare("-per-chromo") == 0 ||
	  lower_version.compare("-perchromo") == 0  ||
	  lower_version.compare("--perchromo") == 0){

	  per_chromo = 1;


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

      }else{
	  printf("Argument %s is not a valid argument.\n", args[i].c_str());
	  print_help(interp);
          return TCL_ERROR;
      }

  }
  struct pio_file_t input;


	if(max_per_file != 0 && per_chromo){
		RESULT_LIT("Cannot use max per file option and per chromosome option together.\n");
		return TCL_ERROR;
	}



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
  read_and_write_to_csv(&input, output_basename.c_str(),
  			 col_size, max_per_file,  
  			 (int)bin_format, per_chromo);
  			 
  return TCL_OK;




}
