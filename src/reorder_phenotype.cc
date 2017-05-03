#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <string>
#include <iterator>
#include <map>
#include "solar.h"
using namespace std;
static bool return_isspace (char  x){
	if(isspace(x)){
		return true;
	}
	return false;
}
static vector<string> parse_line(string line){

	line.erase(remove_if(line.begin(), line.end(), return_isspace), line.end());

	vector<string> parsed_line;

	string current_string;
	char current_char;

	for(string::iterator it = line.begin(); it != line.end(); it++){
		current_char = *it;
		if(current_char == ','){
			parsed_line.push_back(current_string);
			current_string.clear();
		}else{
			current_string += current_char;
		}
	}

	parsed_line.push_back(current_string);

	return parsed_line;


}

static vector<string> get_id_list(ifstream & pedindex_in){
	string line;
	getline(pedindex_in, line);
	vector<string> title_line = parse_line(line);
	vector<string> id_list;
	int id_index = -1;
	for(vector<string>::iterator it = title_line.begin() ; it != title_line.end() ; it++){
		string current_string = *it;
		transform(current_string.begin(), current_string.end(),current_string.begin(), ::toupper);
		if(current_string == "ID"){
			id_index = distance(title_line.begin(), it);
			break;
		}
	}

	if(id_index == -1){
		cout << "ID field in pedindex file could not be found.\n";
		return id_list;
	}

	while(getline(pedindex_in,line)){
		id_list.push_back(parse_line(line).at(id_index));
	}

	return id_list;
}

static void reorder_phenotype(string pedindex_file_name, string phenotype_file_name, string output_file_name, string header_filename){
	vector<string> reordered_lines;
	vector<string> id_list;
	bool is_header = false;
	vector<string> header_list;
	if(header_filename != ""){
		 is_header = true;
		 
		 ifstream header_in(header_filename.c_str());
		 if(header_in.is_open() == false){
			 cout << "Error in opening header file: " << header_filename << "\n";
			 return;
		 }
		 
		 string current_trait;
		 while (header_in >> current_trait){
			 header_list.push_back(current_trait);
		 }
		 
		 header_in.close();
		 
	 }
		
	ifstream pedindex_in(pedindex_file_name.c_str());

	if(pedindex_in.is_open() == false){
		cout << "Error in opening pedindex file: " << pedindex_file_name << "\n";
		return;
	}

	id_list = get_id_list(pedindex_in);
	pedindex_in.close();
	if(id_list.size() == 0)
		return;


	ifstream phenotype_in(phenotype_file_name.c_str());

	if(phenotype_in.is_open() == false){
		cout << "Error in opening phenotype file: " << phenotype_file_name << "\n";
		pedindex_in.close();
		return ;
	}
	map<string, string > id_map;
	string line;
	getline(phenotype_in,line);
	vector<string> title_line = parse_line(line);
	vector<int> index_list;
	int id_index = -1;
	
	ofstream header_out;
	ofstream new_header;
	string header_out_filename =  output_file_name + ".mat.csv";
	string new_header_out_filename = output_file_name + ".header";
	if(is_header){
		header_out.open(header_out_filename.c_str(), ofstream::out);
		new_header.open(new_header_out_filename.c_str(), ofstream::out);
	}
	
	for(vector<string>::iterator it = title_line.begin() ; it != title_line.end(); it++){
		vector<string>::iterator header_it = find(header_list.begin(), header_list.end(), *it);
		if(header_it != header_list.end()){
			index_list.push_back(distance(title_line.begin(), it)); 
			new_header << *header_it << " ";
		}
		string current_string = *it;
		transform(current_string.begin(), current_string.end(),current_string.begin(), ::toupper);
		if(current_string == "ID"){
			id_index = distance(title_line.begin(), it);
			if(is_header == false)
				break;
		}
	}
	if(is_header)
		new_header.close();

	if(id_index == -1){
		cout << "Error in searching for id field in phenotype file:" << phenotype_file_name << "\n";
		phenotype_in.close();
		return;
	}

	ofstream output(output_file_name.c_str());

	output << line << "\n";
	

	map<string, vector<string> > header_map;
	while(getline(phenotype_in, line)){
		vector<string> parsed_line = parse_line(line);
		id_map.insert(pair<string, string >(parsed_line[id_index], line));
		
		if(is_header){
			vector<string> line_out;
			
			for(vector<int>::iterator it_2 = index_list.begin(); it_2 != index_list.end(); it_2++){
				if((it_2 + 1) == index_list.end()){
					string out_value = parsed_line[*it_2];
					if(out_value == ""){
						line_out.clear();
						break;
					}
					line_out.push_back(parsed_line[*it_2] + "\n");
					
				}else{
					string out_value = parsed_line[*it_2];
					if(out_value == ""){
						line_out.clear();
						break;
					}
					line_out.push_back(parsed_line[*it_2] + ",");
				}
			}
			if(line_out.size() != 0){
				header_map.insert(pair<string, vector<string> >(parsed_line[id_index], line_out));
			}
			
		}		
	}
		


	phenotype_in.close();

	for(vector<string>::iterator it = id_list.begin(); it != id_list.end(); it++){
		map<string, string>::iterator map_it = id_map.find(*it);
		if(map_it != id_map.end()){
			output << id_map[*it].data() << "\n";
			if(header_map.find(*it) != header_map.end()){
				for(vector<string>::iterator iter = header_map[*it].begin();  iter != header_map[*it].end(); iter++){
					header_out << *iter;
				}
			}
		}
			

	}
	output.close();
	if(is_header)
		header_out.close();
	return;
}

void print_help(){
	cout << "reorder_phenotype\n";
	cout << "Reorders the ids of a phenotype so that the correspond to that of the pedindex file.\n";
	cout << "reorder_phenotype -ped <pedindex.csv file> -pheno <phenotype file name> -output <output file name> optional: -header <header file name>\n";
}

extern "C" int  RunreorderphenoCmd(ClientData clientData, Tcl_Interp *interp, 
						int argc , const char * argv[]){
	string ped_filename;
	string pheno_filename;
	string output_filename;
	string header_filename;
	for (int index = 1 ; index < argc ; index++){
		string current_arg(argv[index]);
		string upper_arg = current_arg;
		transform(upper_arg.begin(), upper_arg.end(),upper_arg.begin(), ::toupper);
		if((upper_arg == "--HELP") || (upper_arg == "--H") || (upper_arg == "-HELP") || (upper_arg == "-H")){
			print_help();
			return TCL_OK;
		}
	}

	for (int index = 1 ; index < argc ; index += 2){
		string current_arg(argv[index]);
		string upper_arg = current_arg;
		transform(upper_arg.begin(), upper_arg.end(),upper_arg.begin(), ::toupper);
		if(((upper_arg == "--PED") || (upper_arg == "-PED") || (upper_arg == "--PEDIGREE") || (upper_arg == "-PEDIGREE")) && (index + 1 < argc)){
			ped_filename = argv[index + 1];
		}else if (((upper_arg == "--PHENO") || (upper_arg == "-PHENO") || (upper_arg == "--PHENOTYPE") || (upper_arg == "-PHENOTYPE")) && (index + 1 < argc)){
			pheno_filename = argv[index + 1];
		}else if (((upper_arg == "--OUT") || (upper_arg == "-OUT") || (upper_arg == "--OUTPUT") || (upper_arg == "-OUTPUT") ||
				(upper_arg == "-O") || (upper_arg == "--O")) && (index + 1 < argc)){
			output_filename = argv[index + 1];
		}else if (((upper_arg == "--HEAD") || (upper_arg == "-HEAD") || (upper_arg == "--HEADER") || (upper_arg == "-HEADER")) && (index + 1 < argc)){
			header_filename = argv[index + 1];
		}else{
			cout << "Error in argument: "<< argv[index] << "\n";
			print_help();
			return TCL_ERROR;
		}
	}
	
	


	reorder_phenotype(ped_filename,  pheno_filename,  output_filename, header_filename);


	return TCL_OK;
}






