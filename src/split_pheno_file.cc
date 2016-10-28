#include "solar.h"
#include <stdio.h>
#include <stdlib.h>
#include "safelib.h"
#include  <vector>
#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <sstream>

extern std::vector<std::string> parse_line(std::string currentLine);
static bool return_isspace (char  x){
	if(std::isspace(x)){
		return true;
	}
	
	
	return false;
	
}

static void print_vector_line(std::vector<std::string> current_vector, std::ofstream& current_stream){
	
	for(size_t i = 0 ; i < current_vector.size(); i++){
		
		if(i != current_vector.size() - 1){
			current_stream << current_vector[i] << ",";
		}else{
		   current_stream << current_vector[i] << "\n";
	   }
	}
	
}

static void print_vector_line_2(std::vector<std::string> current_vector, std::ostream& current_stream){
	
	for(size_t i = 0 ; i < current_vector.size(); i++){
		
		if(i != current_vector.size() - 1){
			current_stream << current_vector[i] << ",";
		}else{
		   current_stream << current_vector[i] << "\n";
	   }
	}
	
}
std::vector<std::string> parse_line_2(std::string currentLine){
	
	currentLine.erase(std::remove_if(currentLine.begin(), currentLine.end(), return_isspace), currentLine.end());
	
	std::vector<std::string> output_vector;
	std::stringstream ss(currentLine);
	char c;
	std::string curr_str = "";
	while (ss >> c){
		if(c == ','){
			output_vector.push_back(curr_str);
			curr_str = "";
		}else{
			curr_str += c;
		}
		
	}
	output_vector.push_back(curr_str);
	
	return output_vector;
		 
} 



extern "C" int RunSplitPhenoFileCmd (ClientData clientData, Tcl_Interp* interp,
			  int argc, char* argv[]) {
				  
	std::vector<std::string> arglist;
	for(int i = 0; i < argc ; i++){
		std::string currArg(argv[i]);
		arglist.push_back(currArg);
	}
	
	
	std::ifstream file_in(argv[1]);
	
	if(!file_in.is_open()){
		printf("Error: File %s not found\n", argv[1]);
		return TCL_ERROR;
	}
	
	std::vector < std::vector<std::string> >filelines;
	
	std::string eachLine;
	std::getline(file_in, eachLine);
	filelines.push_back(parse_line_2(eachLine));
	eachLine.clear();
    while(std::getline(file_in, eachLine)){

		filelines.push_back(parse_line_2(eachLine));
		
		eachLine.clear();
	}	
		

	std::vector<std::string> titleline(filelines[0]);
	filelines.erase(filelines.begin());
	int class_index = -1;
	for(int i = 0 ; i < titleline.size(); i++){
		
		std::string curr_str = titleline[i];
		if(!StringCmp(curr_str.c_str(), "class", case_ins)){
			class_index = i;
			break;
		}
		
	}
	
	
	if(class_index == -1){
		printf("Error: Class column not found in file\n");
		return TCL_ERROR;
	}
	
	titleline.erase(titleline.begin() + class_index);
	while(filelines.size() != 0){
		std::vector<std::string> front_line = filelines[0];
		
		filelines.erase(filelines.begin());
		std::string current_class;

		current_class = front_line[class_index];
		
		if(class_index >= front_line.size()){
			printf("Error: Column missing in file\n");
			return TCL_ERROR;
		}
		front_line.erase(front_line.begin() + class_index);
		std::string filename(argv[1]);
		filename  = current_class + '_' + filename;
		std::ofstream out_file(filename.c_str());
		print_vector_line(titleline, out_file);
		
		print_vector_line(front_line, out_file);
		
		size_t iterator = 0;
		
		while(iterator < filelines.size()){

			std::vector<std::string> current_line(filelines[iterator]);
			if(class_index >= current_line.size()){
				printf("Error: Column missing in file\n");
				return TCL_OK;
			}
			
			if(!StringCmp(current_line[class_index].c_str(), current_class.c_str(), case_sens)){
				
				current_line.erase(current_line.begin() + class_index);
				print_vector_line(current_line, out_file);
				filelines.erase(filelines.begin() + iterator);
				
			}else{
				iterator += 1;
			}
			
		}
		
		out_file.close();
		
	}
	
	
	
	return TCL_OK;
				
	
}
	
	
