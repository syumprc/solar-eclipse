#include "RicVolumeSet.h"
#include "solar.h"
#include <string>
#include <vector>
#include <queue>
using namespace std;



extern "C" void cdfnor_ (int* which, double* p, double* q, double* x, double* mean, double* sd, int* status, double* bound);

static inline double compute_inverse_normal(double  pct){
	
	double q = 1.0 - pct;
	double mean = 0.0;
	double standard_deviation = 1.0;
	int status = 0;
	int which = 2;
	double bound = 0;
	double x = 0;
	cdfnor_(&which, &pct, &q, &x, &mean, &standard_deviation, &status, &bound);
	return x;
}

static  void inorm(queue<double> sorted_input_data, queue<string> sorted_input_ids, vector<double> output_data, vector<string> output_ids){
	double current, last = -1.0;
	double pct, z;
	const size_t count = input_data.size();
	size_t position = 0, shared_count;
	string current_id;
	vector<double> copies;
	vector<string> copy_ids;
	double shared_value, sum = 0.0;
	while(input_data.size() != 0 ){
		shared_value = sorted_input_data.front();
		while(shared_value == sorted_input_data.front()){
			current = sorted_input_data.pop();
			++shared_count;
			pct = double(position++)/(count + 1);
			z = compute_inverse_normal(pct);
			sum += z;
		}
		sum /= shared_count;
		for(size_t id = 0; id < shared_count; id++){
			current_id = sort_input_ids.pop();
			output_data.push_back(sum);
			output_ids.push_back(current_id);
		}
		
		sum = 0.0;
		shared_count  = 0;
		
		
	}
}

static void 

