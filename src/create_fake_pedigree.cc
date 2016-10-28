#include "solar.h"
#include <stdio.h>
#include <stdlib.h>
#include "safelib.h"
#include  <vector>
#include <fstream>
#include <string>
#include <string.h>
#include <iterator>
#include <algorithm>
#include <iostream>



void printHelp(){
	printf("Purpose: This command creates a pedigree file given a phenotype file taken as input.\n \n");
	printf("Usage: create_fake_pedigree <phenotype filename> [-o output pedigree filename\n");
	printf("  <phenotype filename> Phenotype filename to be used to create pedigree\n");
	printf("  [-o <output pedigree filename>] Option to name output pedigree filename\n");

}




extern std::vector<std::string> parse_line_2(std::string currentLine);

std::vector<std::string> parse_line(std::string currentLine){
	
	std::vector<std::string> output_vector;
	int commas = 0;
	char lastchar = 'a';
	std::string curr_str;
	for(int i = 0; i < currentLine.size(); i++){
		char curr_char = currentLine[i];
		if(curr_char == '\n' || curr_char == '\r'){
			output_vector.push_back(curr_str);
			curr_str.clear();
			
			
			
			continue;
		}
		if(curr_char == ','){
			if(((lastchar == ' ') || (lastchar == ',')) && curr_char ==','){
				output_vector.push_back("");
			}else{
				output_vector.push_back(curr_str);
				curr_str.clear();
			}
		}else{
			curr_str += curr_char;
		}
		lastchar = curr_char;
		
	}
	
	return output_vector;
	
	
} 


std::string convert_line(std::string line_in, bool is_first_line){
	std::vector<std::string> divided_line = parse_line_2(line_in);
	
	if(is_first_line){
		divided_line.push_back("fa");
		divided_line.push_back("mo");
	}else{
		divided_line.push_back("0");
		divided_line.push_back("0");
	}
	
	std::string output_string = "";
	
	for(int i = 0; i < divided_line.size(); i++){
		if(i == divided_line.size() - 1){
			output_string += divided_line[i];
		}else{
			output_string += divided_line[i] + ",";
		}
		
	}
	
	return output_string;
	
	
}



extern "C" int RunCreateFakePedigreeCmd(ClientData clientData, Tcl_Interp *interp,
int argc, char *argv[]){
	
	if(argc != 2 && argc != 4){
		printf("Error: Arguments are incorrect see help\n");
		return TCL_ERROR;
	}
	
	if(!StringCmp(argv[1], "help", case_ins)){
		printHelp();
		return TCL_OK;
	}
	std::string output_name(argv[1]);
	if(argc == 4){
		if(!StringCmp(argv[2], "-o", case_ins)){
			std::string temp(argv[3]);
			output_name = temp;
		}else {
			printf("Error: Arguments are incorrect see help\n");
			return TCL_ERROR;
		}
	}else{
		output_name = "pedigree_" + output_name;
	}
	
	std::ifstream file_in(argv[1]);
	if(!file_in.is_open()){
		printf("Error: File %s not found\n", argv[1]);
		return TCL_ERROR;
	}
	
	
	
	
	std::vector<std::string> filelines;
	
	std::string eachLine;
	std::getline(file_in, eachLine);
	filelines.push_back(convert_line(eachLine, true));
	eachLine.clear();
    while(std::getline(file_in, eachLine)){

		filelines.push_back(convert_line(eachLine, false));
		
		eachLine.clear();
	}	
	
	
	file_in.close();
	
	
	
	std::ofstream file_out(output_name.c_str());
	for(int i = 0 ; i < filelines.size(); i++){
		file_out << filelines[i] << "\n";
	}
	
	file_out.close();
	std::string run_command = "load ped " + output_name; 
	int success = Solar_Eval(interp, run_command.c_str());
	if(success == TCL_ERROR){
		printf("Error in loading newly created pedigree file\n");
		return TCL_ERROR;
	}
	
	
	return TCL_OK;
	
}

	
	
	
	
	 
	
	
