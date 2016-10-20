#include<cmath>
#include<stdio.h>
#include<vector>
#include<algorithm>
#include<sstream>
#include<fstream>
#include<iostream>
#include<string>
#include "solar.h"
#include <time.h>
static bool return_isspace (char  x){
	if(std::isspace(x)){
		return true;
	}
	return false;
}

static std::vector<float> parse_line(std::string currentLine){

	currentLine.erase(std::remove_if(currentLine.begin(), currentLine.end(), return_isspace), currentLine.end());
	std::vector<float> output_vector;
	std::stringstream ss(currentLine);
	char c;
	std::string curr_str = "";
	while (ss >> c){
		if(c == ','){
		    output_vector.push_back(std::atof(curr_str.c_str()));
		    curr_str = "";
		}else{
			curr_str += c;
		}

	}
	output_vector.push_back(std::atof(curr_str.c_str()));

	return output_vector;

}


/*
 * Input: filename (file to be opened), rows (pointer to the number of rows in filename),
 * and columns (number of columns in filename)
 * Output: A float vector containing all the values in filename (vector is empty if there
 * is an error).
 */
static float * csvread_float(const char * filename, size_t * rows,
	    size_t * columns, bool column_major) {

	float * data;
	std::ifstream inputData(filename);
	if(inputData.is_open() == false){
		std::cerr << "Error in opening file: " << filename << "\n";
		return data;
	}
	std::string Line;

	std::vector< std::vector<float> >   lines;
	while (std::getline(inputData, Line)) {
		lines.push_back(parse_line(Line));
	}
	inputData.close();
	*rows = lines.size();
	*columns = lines[0].size();
	data = new float[*rows*(*columns)];
	if(column_major == true){
	    for(size_t row = 0 ; row < *rows ; row++){
		std::vector<float> temp(lines[row]);
		for(size_t column = 0 ; column < *columns ; column++){
		    data[column*(*rows) + row] = temp[column];
		}
	    }
	}else{
	    for(size_t row = 0 ; row < *rows ; row++){
		std::vector<float> temp(lines[row]);
		for(size_t column = 0 ; column < *columns ; column++){
		    data[row*(*columns) + column] = temp[column];
		}
	    }
	}



	return data;

}




static int random_num( int numMax, int numMin) {

	//float num = (float) rand() / RAND_MAX;
	return numMin + rand()%numMax;
}


static unsigned int * getPvector(unsigned int n_subjects, unsigned int n_permutations) {

  unsigned int * Pvec = new unsigned int[n_subjects * n_permutations];
	for(unsigned int j = 0 ; j < n_permutations; j++){
		std::vector<unsigned int> randNums;

		for (int r = 0; r < n_subjects; r++) {
			randNums.push_back(r);
		}

	    for (unsigned int n = n_subjects; n > 0; n--) {

		unsigned int Idx = rand()%n;
		Pvec[(n - 1) + j*n_subjects] = randNums[Idx];
		randNums.erase(randNums.begin() + Idx);
	    }
	   // Pvec[j] = randNums[0];

	}






	return Pvec;

}

static std::string convert_float_to_string(float value){
  std::stringstream ss;
  ss << value;
  std::string str_value(ss.str());

  return str_value;
}
static void write_to_file(std::string file_name, float * h2, float * indicator, float * pval, bool get_pvalue, size_t n_voxels){

  std::ifstream file_in(file_name.c_str());
  std::string * lines = new std::string[n_voxels];


  for(size_t voxel = 0; voxel < n_voxels; voxel++){
      
      std::getline(file_in,lines[voxel]);

      lines[voxel] += ',' +  convert_float_to_string(indicator[voxel]) + ',' + convert_float_to_string(h2[voxel]);
      if(get_pvalue){
		  lines[voxel] +=  ',' + convert_float_to_string(pval[voxel]) + '\n';
      }else{
		  lines[voxel] += '\n';
      }
  }

  file_in.close();

  std::ofstream file_out(file_name.c_str());
  file_out << "Trait,Indicator,H2r";
  if (get_pvalue){
	  file_out << ",pvalue\n";
  }else{
	  file_out << "\n";
  }
  for(size_t voxel = 0; voxel < n_voxels; voxel++){
      file_out << lines[voxel].c_str();
  }

  file_out.close();
  delete [] lines;


}
extern "C" std::vector<int> select_devices();
extern "C" int call_cudafPHI(float * h2, float * indicator, float * pvals,float * h_y,float * h_Z,float * h_hat,
    float * h_evectors, unsigned int * h_pmatrix, bool covariates, bool get_pval,
    size_t n_voxels, size_t n_subjects, size_t n_permutations, std::vector<int> selected_devices);


extern "C" int RuncudafphiCmd(ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[]){
			  
			  
  std::vector<int> selected_devices = select_devices();

  if(selected_devices.size() == 0){
      printf("No usable devices were selected.\n");
      return EXIT_FAILURE;
  }
			  
  size_t  n_subjects;

  size_t  n_permutations;
  size_t  n_voxels;
  float * h_y = csvread_float(argv[1], &n_subjects, &n_voxels, true);
  float * e_vectors = csvread_float(argv[3], &n_subjects, &n_subjects, true);

  srand(time(NULL));

  size_t  dummy;
  float * hat;
  float *  aux = csvread_float(argv[2], &n_subjects, &dummy, true);

  n_permutations = std::atoi(argv[5]);



  unsigned int * p_matrix;
  bool covariates = false;
  int get_pval_int;
  if(argc == 8){
      hat = csvread_float(argv[argc - 1], &n_subjects, &n_subjects, false);
      covariates = true;
      get_pval_int = std::atoi(argv[argc - 2]);
   }else{
       get_pval_int = std::atoi(argv[argc - 1]);

   }
  srand(time(NULL));
  bool get_pval = false;
  if(get_pval_int == 1){
      get_pval = true;
      p_matrix = getPvector(n_subjects, n_permutations);
  }

  std::string output_filename(argv[4]);

  float * indicator = new float[n_voxels];
  float * pvals = new float[n_voxels];
  float * h2 = new float[n_voxels];
  call_cudafPHI(h2, indicator, pvals, h_y, aux, hat, e_vectors, p_matrix, covariates, get_pval, n_voxels, n_subjects, n_permutations, selected_devices);

  write_to_file(output_filename, h2, indicator, pvals, get_pval, n_voxels);
  delete [] h2;
  delete [] pvals;
  delete [] indicator;
  delete [] aux;
  delete [] e_vectors;
  delete [] h_y;
  if(covariates)
    delete [] hat;
  if(get_pval)
    delete [] p_matrix;


  return EXIT_SUCCESS;
}
