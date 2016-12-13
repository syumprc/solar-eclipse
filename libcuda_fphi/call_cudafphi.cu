#include "cudafphi.cuh"
#include <iostream>
#include <cuda_runtime.h>
#include<stdio.h>
#include<time.h>
#include<pthread.h>
#include<vector>
#include <algorithm>
#include <iterator>
#include <chrono>
#include <Eigen/Dense>
#include <fstream>
#include <curand.h>
#include <curand_kernel.h>
#include <thread>

#define pthreadErrchk(ans) { pthreadAssert((ans), __FILE__, __LINE__); }
static inline void pthreadAssert(int code, const char *file, int line, bool abort=true)
{

   if (code != cudaSuccess)
   {
	   fprintf(stderr,"pthreadAssert: %d %s %d\n", code, file, line);

	   fprintf(stdout,"pthreadAssert: %d %s %d\n", code, file, line);

	   if (abort) exit(code);
   }else{
	   fprintf(stderr,"pthreadAssert: %d %s %d\n", code, file, line);
   }
}



/*
void inverse(float* A, int N)
{
    int *IPIV = new int[N+1];
    int LWORK = N*N;
    float *WORK = new float[LWORK];
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    delete IPIV;
    delete WORK;
}


*/
size_t iterations;
size_t batch_size;
size_t free_voxels;

size_t total_voxels;


static __global__ void Inv4by4(float * ZTZI) {
	float a = ZTZI[0];
	float b = ZTZI[1];
	float c = ZTZI[2];
	float d = ZTZI[3];

	float det = a * d - b * c;
	ZTZI[0] = d / det;
	ZTZI[1] = -c / det;
	ZTZI[2] = -b / det;
	ZTZI[3] = a / det;

}
static __global__ void  set_Identity(float * I, size_t n_subjects){
	size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;
	size_t colIdx = threadIdx.y + blockIdx.y*blockDim.y;

	if(rowIdx < n_subjects && colIdx < n_subjects){
		if(colIdx == rowIdx){
			I[colIdx*n_subjects + rowIdx] = 1.f;
		}else{
			I[colIdx*n_subjects + rowIdx] = 0.f;
		}
	}

}
/*
static __global__ void calculate_U(float * d_A, int i, size_t n_covariates){

	int k = threadIdx.x+ i + 1;
	int j = threadIdx.y + i + 1;
	const float value = d_A[k*n_covariates + j];

	d_A[k*n_covariates + j] = value - d_A[i + k*n_covariates]*d_A[j + i*n_covariates];





}

static __global__ void calculate_L(float * d_A, int i, size_t n_covariates){


	int j = threadIdx.x + 1;
	const float value = d_A[i*n_covariates + j];
	d_A[i*n_covariates + j] = value/d_A[i*n_covariates + i];

	__syncthreads();

	if(j == 0){
		dim3 blockSize(n_covariates - i - 1, n_covariates - i - 1);
		calculate_U<<<1, blockSize, 0, 0>>>(d_A, i, n_covariates);
	}



}

static __global__ void solve_L(const float * d_A, float * d_B, size_t n_covariates){



}
*/

static __global__ void init_rand_kernel(curandState * state, unsigned int seed, size_t n_permutations){

	const size_t pIdx = threadIdx.x + blockDim.x*blockIdx.x;

	if(pIdx < n_permutations){
		curand_init(seed, pIdx, 0, &state[pIdx]);
	}

}

static __device__ unsigned int generate_number(curandState * state, unsigned int  max){

	return curand(state) % max;

}

static __global__ void init_permutation_matrix(unsigned int * pmatrix, unsigned int * max, size_t n_subjects, size_t n_permutations){

	const size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;
	const size_t pIdx = threadIdx.y + blockDim.y*blockIdx.y;

	if(rowIdx >= n_subjects || pIdx >= n_permutations)
		return;



	pmatrix[rowIdx + pIdx*n_subjects] = n_subjects + 1;

	if(rowIdx == 0)
		max[pIdx] = 1;


}

static __global__ void create_permutation_matrix(curandState * state, unsigned int * pmatrix, unsigned int * max, size_t n_subjects , size_t n_permutations){


	const size_t rowIdx = n_subjects - 1 - (threadIdx.x + blockIdx.x*blockDim.x);
	const size_t pIdx = threadIdx.y + blockDim.y*blockIdx.y;

	if(rowIdx >= n_subjects || pIdx >= n_permutations)
		return;



	curandState  l_state = state[pIdx];




	unsigned int Idx;
	unsigned int current_Idx;

	do{
		Idx = generate_number(&l_state, max[pIdx]);
		state[pIdx] = l_state;
		current_Idx = atomicCAS(&pmatrix[Idx + pIdx*n_subjects], n_subjects + 1, Idx);

		if(current_Idx != n_subjects  + 1){
			pmatrix[rowIdx + pIdx*n_subjects] = Idx;
			atomicInc(&max[pIdx], 1);
		}
		__threadfence();


	}while(current_Idx != n_subjects + 1);
}


static __global__ void  compute_Inverse(float * d_A, size_t n_covariates){

	  float * h_A = &d_A[0];

	  for(size_t i = 0 ; i < n_covariates ; i++){
	      for(size_t j = 0 ; j < i; j++){
	    	  float sum = h_A[j*n_covariates + i];
	    	  for(size_t k = 0 ; k < j ;k++){
	    		  sum -= h_A[k*n_covariates + j]*h_A[j*n_covariates + k];
	    	  }
	    	  h_A[j*n_covariates + i] = sum/h_A[j*n_covariates + j];
	      }
	      for(size_t j = i ; j < n_covariates; j++){
	    	  float sum = h_A[j*n_covariates + i];
	    	  for(size_t k = 0 ; k < i ;k++){
	    		  sum -= h_A[k*n_covariates + j]*h_A[j*n_covariates + k];
	    	  }
	    	  h_A[j*n_covariates + i] = sum;
	      }

	  }


	  float * M = new float[n_covariates*n_covariates];

	  for(size_t j = 0 ; j < n_covariates ; j++){
	      for(size_t i = 0 ; i < n_covariates ; i++){
		  float sum = 0.0;
		  if(i == j)
		    sum = 1.0;
		  for(size_t k = 0 ; k < i; k++){
		      sum -= M[j*n_covariates + k]*h_A[k*n_covariates + i];
		  }
		  M[j*n_covariates + i] = sum;
	      }
	  }




	  float * h_AI = new float[n_covariates*n_covariates];
	  for(size_t j = 0;  j < n_covariates; j++){
	      for(size_t i = 0; i < n_covariates ; i++){
		  float sum = M[j*n_covariates + i];
		  for(size_t k = 0; k < i ; k++){
		      sum -= h_AI[j*n_covariates + k]*h_A[k*n_covariates + i];
		  }
		  h_AI[j*n_covariates + i] = sum/h_A[i*n_covariates + i];
	      }
	  }

	delete [] M;


	for(int j = 0 ; j < n_covariates; j++){
		for(int k = 0 ; k < n_covariates ; k++){
			d_A[j*n_covariates +k] = h_AI[j*n_covariates + k];
		}
	}
	delete [] h_AI;
}



static __global__ void  LU_factor(float * L, float * U,  float * d_A, size_t n_covariates){




	for(int k = 0 ; k < n_covariates ; k++){
		L[k*n_covariates + k] = 1.f;
		for(int j = k ; j < n_covariates; j++){
			 float sum = 0.f;

			for(int s = 0 ; s < k - 1; s++){
				 float L_val = L[s*n_covariates + k];
				 float U_val = U[j*n_covariates + s];
				sum += L_val*U_val;
			}

			U[j*n_covariates + k] = d_A[j*n_covariates + k] - sum;

		}

		for(int i = k+1; i < n_covariates ; i++){
			 float sum = 0;

			for(int s = 0 ;s < k  - 1 ;s++){
				 float L_val = L[s*n_covariates + i];
				 float U_val = U[k*n_covariates + s];
				sum += L_val*U_val;
			}

			L[k*n_covariates + i] = (d_A[k*n_covariates + i] - sum)/U[k*n_covariates + k];
		}
	}



}

static __global__ void subtract_matrices(float * d_hat, size_t n_subjects){

	size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;

	size_t colIdx = threadIdx.y + blockIdx.y*blockDim.y;

	if(rowIdx < n_subjects && colIdx < n_subjects){
		if(rowIdx == colIdx){
			d_hat[colIdx*n_subjects + rowIdx] = 1.f - d_hat[colIdx*n_subjects + rowIdx];
		}else{
			d_hat[colIdx*n_subjects + rowIdx] = -d_hat[colIdx*n_subjects + rowIdx];
		}
	}


}

static void compute_hat_CPU(float * h_cov, float * h_hat, size_t n_subjects, size_t n_covariates){
	Eigen::MatrixXf cov(n_subjects, n_covariates);

	for(size_t row = 0 ; row < n_subjects ; row++){
		for(size_t col = 0 ; col < n_covariates ;col++){
			cov(row, col)  = h_cov[col*n_subjects + row];
		}
	}
	Eigen::MatrixXf identity(n_covariates,n_covariates);
	for(size_t row = 0 ; row < n_covariates ; row++){
		for(size_t col = 0 ; col < n_covariates ;col++){
			if(col == row){
				identity(row, col) = 1.f;
			}else{
				identity(row, col) = 0.f;
			}
		}
	}
	Eigen::MatrixXf inverse_cov = ((cov.transpose())*cov).ldlt().solve(identity);

	Eigen::MatrixXf big_identity(n_subjects,n_subjects);
	for(size_t row = 0 ; row < n_subjects ; row++){
		for(size_t col = 0 ; col < n_subjects ;col++){
			if(col == row){
				big_identity(row, col) = 1.0;
			}else{
				big_identity(row,col) = 0.f;
			}
		}
	}


	Eigen::MatrixXf intermediate_hat = cov*inverse_cov;


	Eigen::MatrixXf hat = (intermediate_hat*cov.transpose());


	hat = big_identity - hat;
	for(size_t row = 0 ; row < n_subjects ; row++){
		for(size_t col = 0 ; col < n_subjects ;col++){
			h_hat[col*n_subjects + row] = hat(row, col);
		}
	}


}

static void compute_new_hat_CPU(float * h_cov, float * h_hat , float * h_sy, size_t n_subjects , size_t n_covariates, size_t cov_number){
	Eigen::MatrixXf cov(n_subjects, n_covariates + 1);

	for(size_t row = 0 ; row < n_subjects ; row++){
		for(size_t col = 0 ; col < n_covariates + 1 ;col++){
			if(col == (n_covariates)){
				cov(row, col)  = h_sy[cov_number*n_subjects + row];
			}else{
				cov(row, col)  = h_cov[col*n_subjects + row];
			}
		}
	}
	Eigen::MatrixXf identity(n_covariates + 1,n_covariates+1);
	for(size_t row = 0 ; row < n_covariates + 1 ; row++){
		for(size_t col = 0 ; col < n_covariates + 1 ;col++){
			if(col == row){
				identity(row, col) = 1.0;
			}else{
				identity(row, col) = 0.f;
			}
		}
	}
	Eigen::MatrixXf inverse_cov = ((cov.transpose())*cov).ldlt().solve(identity);

	Eigen::MatrixXf big_identity(n_subjects,n_subjects);
	for(size_t row = 0 ; row < n_subjects ; row++){
		for(size_t col = 0 ; col < n_subjects ;col++){
			if(col == row){
				big_identity(row, col) = 1.0;
			}else{
				big_identity(row,col) = 0.f;
			}
		}
	}


	Eigen::MatrixXf intermediate_hat = cov*inverse_cov;


	Eigen::MatrixXf hat = (intermediate_hat*cov.transpose());


	hat = big_identity - hat;

	for(size_t row = 0 ; row < n_subjects ; row++){
		for(size_t col = 0 ; col < n_subjects ;col++){
			h_hat[col*n_subjects + row] = hat(row, col);
		}
	}


}

static void compute_ZTZI(const float * h_Z, float * h_ZTZI, size_t n_subjects){

	Eigen::MatrixXf Z(n_subjects, 2);

	for(size_t col = 0 ; col < 2 ; col++){
		for(size_t row = 0 ; row < n_subjects ; row++){
			Z(row, col) = h_Z[col*n_subjects + row];
		}
	}
	Eigen::MatrixXf ZTZI = (Z.transpose()*Z);

	 float det = ZTZI(0, 0)*ZTZI(1,1) - ZTZI(1, 0)*ZTZI(0, 1);
	h_ZTZI[0] =  ZTZI(1, 1)/det;
	h_ZTZI[1] = -ZTZI(0, 1)/det;
	h_ZTZI[2] = -ZTZI(1, 0)/det;
	h_ZTZI[3] = ZTZI(0, 0)/det;



}

extern "C" int call_cudafphi_qsub(float * conn_h2, float * conn_indicator,
		 float * conn_pval, float * h_y,  const float * h_Z,  float * h_cov,
            const float * h_evectors,const unsigned int * h_pmatrix, bool covariates, bool get_pval, bool get_connectivity, int node_number,
            size_t n_voxels, size_t n_subjects, size_t n_permutations, size_t qsub_batch_size, size_t qsub_shared_batch_size,
            size_t n_covariates, std::vector<int> selected_devices){

 int devices = selected_devices.size();
 if(n_voxels <= BLOCK_SIZE_3){
   batch_size = n_voxels;
   iterations = 1;
   free_voxels = batch_size;

 }else{
   batch_size = BLOCK_SIZE_3;
   iterations = ceil(float(n_voxels)/float(batch_size));
   free_voxels = n_voxels % batch_size;
   if(n_voxels % batch_size == 0)
     free_voxels = batch_size;
 }

 float * h_Y[devices];



 int n_streams = iterations;
 cuda_fphi_variables_per_stream stream_vars[n_streams];
 cuda_fphi_variables_per_device device_vars[devices];

 cudaStream_t device_stream[devices];

 cublasHandle_t device_handles[devices];
 std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

	int blockSize_n_subjects = 32;

	if(n_subjects >= 64){
		blockSize_n_subjects = 64;
	}
	if(n_subjects >= 128){
		blockSize_n_subjects = 128;
	}
	if(n_subjects >= 256){
		blockSize_n_subjects = 256;
	}
	if(n_subjects >= 512){
		blockSize_n_subjects = 512;
	}
	if(n_subjects >= 1024){
		if(n_subjects % 1024 <= n_subjects % 512){
			blockSize_n_subjects = 1024;
		}else{
			blockSize_n_subjects = 512;
		}

	}

	dim3 blockSize(blockSize_n_subjects, 1, 1);
	dim3 gridSize(ceil(float(n_subjects)/float(blockSize_n_subjects)),n_subjects, 1);
 for(int current_device = 0 ; current_device < devices; current_device++){

	  gpuErrchk(cudaSetDevice(selected_devices[current_device]));
	  device_vars[current_device].device_id = selected_devices[current_device];
	  gpuErrchk(cudaSetDeviceFlags((unsigned int)cudaDeviceScheduleAuto));
	  gpuErrchk(cudaStreamCreate(&device_stream[current_device]));
	  cublasErrchk(cublasCreate_v2(&device_handles[current_device]));
	  cublasErrchk(cublasSetStream_v2(device_handles[current_device], device_stream[current_device]));

 }


 float alpha = 1.f;
 float beta = 0.f;


 float * d_cov[devices];
 float * h_hat = new float[n_subjects*n_subjects];


 for(int current_device = 0; current_device < devices; current_device++){

	  gpuErrchk(cudaSetDevice(device_vars[current_device].device_id));



	  if(covariates){
		  gpuErrchk(cudaMalloc((void**)&device_vars[current_device].hat, sizeof(float)*n_subjects*n_subjects));
	  }

	  gpuErrchk(cudaMalloc((void**)&device_vars[current_device].evectors, sizeof(float)*n_subjects*n_subjects));



     if(get_pval){
    	 gpuErrchk(cudaMalloc((void**)&device_vars[current_device].pmatrix, sizeof(unsigned int)*n_subjects*n_permutations));
     }

     gpuErrchk(cudaMalloc((void**)&device_vars[current_device].aux_vars.d_Z, sizeof(float)*n_subjects*2));



     gpuErrchk(cudaMalloc((void**)&device_vars[current_device].aux_vars.d_ZTZI, sizeof(float)*2*2));



     device_vars[current_device].device_id = selected_devices[current_device];
     device_vars[current_device].get_pval = get_pval;
     device_vars[current_device].covariates = covariates;
     device_vars[current_device].get_connectivity = get_connectivity;

     gpuErrchk(cudaMemcpyAsync((float *)device_vars[current_device].aux_vars.d_Z, h_Z,  sizeof(float)*n_subjects*2, cudaMemcpyHostToDevice, device_stream[current_device]));
     if(get_pval)
    	 gpuErrchk(cudaMemcpyAsync((unsigned int *)device_vars[current_device].pmatrix, h_pmatrix, sizeof(unsigned int)*n_subjects*n_permutations, cudaMemcpyHostToDevice, device_stream[current_device]));
     gpuErrchk(cudaMemcpyAsync((float *)device_vars[current_device].evectors, h_evectors, sizeof(float)*n_subjects*n_subjects, cudaMemcpyHostToDevice, device_stream[current_device]));
     cublasErrchk(cublasSgemm_v2(device_handles[current_device], CUBLAS_OP_T, CUBLAS_OP_N, 2, 2, n_subjects, &alpha, device_vars[current_device].aux_vars.d_Z, n_subjects, device_vars[current_device].aux_vars.d_Z, n_subjects, &beta, device_vars[current_device].aux_vars.d_ZTZI, 2));

     Inv4by4<<<1, 1, 0, device_stream[current_device]>>>(device_vars[current_device].aux_vars.d_ZTZI);

 }

 for(int current_device = 0 ; current_device < devices; current_device++){

	  gpuErrchk(cudaSetDevice(selected_devices[current_device]));

	  gpuErrchk(cudaStreamSynchronize(device_stream[current_device]));

 }

 stream_vars[0].n_permutations = n_permutations;
 stream_vars[0].n_subjects = n_subjects;
 stream_vars[0].shared_device_vars = device_vars[0];

 stream_vars[0].stream_number= 0;
 stream_vars[0].n_voxels = batch_size;

 size_t total_bytes;
 size_t free_bytes;
 unsigned int total_number_of_threads_usable_by_device[devices];
 unsigned int threads_in_use_per_device[n_streams];
 size_t total_usable_threads = 0;
 size_t used_bytes = run_allocation_test(stream_vars[0]);
 for(int device = 0; device < devices ; device++){
     gpuErrchk(cudaSetDevice(selected_devices[device]));
     gpuErrchk(cudaMemGetInfo(&free_bytes, &total_bytes));

     total_number_of_threads_usable_by_device[device] = floor(float(total_bytes)/float(used_bytes));
     total_usable_threads += total_number_of_threads_usable_by_device[device];
     threads_in_use_per_device[device] = 0;

 }

 int device = 0;
 std::vector<int> current_running_streams;
 pthread_t pthreads[n_streams];
 void * return_val;





 float * d_sy;
 float * d_y;
 float * mean;
 gpuErrchk(cudaSetDevice(selected_devices[0]));
 gpuErrchk(cudaMalloc((void**)&d_sy, sizeof(float)*n_subjects*batch_size));
 gpuErrchk(cudaMalloc((void**)&d_y, sizeof(float)*n_subjects*batch_size));
 gpuErrchk(cudaMalloc((void**)&mean, sizeof(float)*batch_size));

 for(size_t stream_number = 0 ; stream_number < iterations ; stream_number++){
	 size_t current_n_voxels;


	 if(stream_number == iterations - 1)
		 current_n_voxels = free_voxels;
	 else
		 current_n_voxels = batch_size;
	 gpuErrchk(cudaMemcpyAsync(d_sy, h_y + stream_number*batch_size*n_subjects, sizeof(float)*n_subjects*current_n_voxels, cudaMemcpyHostToDevice, device_stream[0]));
	 whitened_Y(d_sy, mean,  device_vars[0].evectors,  d_y,  n_subjects,  current_n_voxels, device_handles[0],  device_stream[0]);

	 gpuErrchk(cudaMemcpyAsync(h_y + stream_number*batch_size*n_subjects, d_sy , sizeof(float)*n_subjects*current_n_voxels, cudaMemcpyDeviceToHost, device_stream[0]));

 }
 gpuErrchk(cudaStreamSynchronize(device_stream[0]));
 gpuErrchk(cudaDeviceSynchronize());
 gpuErrchk(cudaFree(mean));
 gpuErrchk(cudaFree(d_y));
 gpuErrchk(cudaFree(d_sy));




 for(size_t current_voxel = 0 ; current_voxel <  qsub_batch_size; current_voxel++){


	 compute_new_hat_CPU(h_cov, h_hat , h_y,  n_subjects ,  n_covariates,  current_voxel + qsub_shared_batch_size*node_number);
	 for(int device_index = 0 ; device_index < devices ; device_index++){
		 gpuErrchk(cudaSetDevice(selected_devices[device_index]));
		 gpuErrchk(cudaMemcpy(device_vars[device_index].hat, h_hat, sizeof(float)*n_subjects*n_subjects, cudaMemcpyHostToDevice));
	 }

	 for(int current_stream = 0 ; current_stream < n_streams; current_stream++){


		 do{
			 if(threads_in_use_per_device[device] >= total_number_of_threads_usable_by_device[device])
				 device++;
			 else
				 break;
		 }while(device != devices);

		 if(device == devices){
			 int streamIdx;
			 do{
				 for(std::vector<int>::iterator it = current_running_streams.begin(); it != current_running_streams.end() ;it++){
					 pthreadErrchk(pthread_tryjoin_np(pthreads[*it], &return_val));
					 if(return_val != NULL){
						 device = std::distance(selected_devices.begin(), std::find( selected_devices.begin(), selected_devices.end(), stream_vars[*it].shared_device_vars.device_id));
						 streamIdx = *it;
						 threads_in_use_per_device[device] -= 1;
						 break;
					 }
	  			  }
			 }while(device == devices);
			 	 current_running_streams.erase( std::find(current_running_streams.begin(), current_running_streams.end(), streamIdx));
			 	 return_val = NULL;
		 }
		 stream_vars[current_stream].shared_device_vars = device_vars[device];
		 stream_vars[current_stream].shared_device_vars.h2 = &conn_h2[current_voxel*n_voxels];
		 if(get_pval)
			 stream_vars[current_stream].shared_device_vars.pvals = &conn_pval[current_voxel*n_voxels];
		 stream_vars[current_stream].shared_device_vars.indicator = &conn_indicator[current_voxel*n_voxels];
		 stream_vars[current_stream].h_y = h_y;
		 stream_vars[current_stream].h_cov = h_cov;

		 stream_vars[current_stream].n_permutations = n_permutations;
		 stream_vars[current_stream].n_subjects = n_subjects;
		 stream_vars[current_stream].n_voxels = batch_size;
		 stream_vars[current_stream].n_covariates = n_covariates;
		 stream_vars[current_stream].stream_number = current_stream;
		 stream_vars[current_stream].cov_number = current_voxel +  qsub_shared_batch_size*node_number;

		 if(current_stream == iterations - 1)
			 stream_vars[current_stream].n_voxels = free_voxels;
		 pthreadErrchk(pthread_create(&pthreads[current_stream], NULL, run_cudafphi_connect_pthread, (void*)&stream_vars[current_stream]));
		 current_running_streams.push_back(current_stream);
		 threads_in_use_per_device[device] += 1;
		 device++;
		 if(device == devices)
			 device = 0;



	 }

	 for(int current_device = 0 ; current_device < devices; current_device++){
		 threads_in_use_per_device[current_device] = 0;
	 }



	 for(int it = 0; it < current_running_streams.size() ; it++){
		 pthreadErrchk(pthread_join(pthreads[current_running_streams[it]], NULL));
	 }



	 current_running_streams.clear();
	 device = 0;

 }



 for(std::vector<int>::iterator device_it = selected_devices.begin() ; device_it != selected_devices.end(); device_it++){
	  gpuErrchk(cudaSetDevice(*device_it));
	  const size_t  device_Idx = std::distance(selected_devices.begin(), device_it);

	  gpuErrchk(cudaStreamDestroy(device_stream[device_Idx]));
	  gpuErrchk(cudaDeviceSynchronize());
	  gpuErrchk(cudaFree((void*)device_vars[device_Idx].evectors));
     if(get_pval){
    	 gpuErrchk(cudaFree((void*)device_vars[device_Idx].pmatrix));
     }
     gpuErrchk(cudaFree((void*)device_vars[device_Idx].aux_vars.d_Z));
     gpuErrchk(cudaFree((void*)device_vars[device_Idx].aux_vars.d_ZTZI));

     cublasErrchk(cublasDestroy_v2(device_handles[device_Idx]));

     gpuErrchk(cudaFree((void*)device_vars[device_Idx].hat));

 }

 delete [] h_hat;


 auto elapsed = std::chrono::high_resolution_clock::now() - start;
 long long seconds = std::chrono::duration_cast<std::chrono::seconds>(elapsed).count();

 printf ("It took %i seconds.\n", (((int)(seconds))));

 return 1;

}

 extern "C" int call_cudafphi(float * h2, float * indicator, float * pvals, float * conn_h2, float * conn_indicator,
		 float * conn_pval, float * h_y,  const float * h_Z,  float * h_cov,
             const float * h_evectors, unsigned int * h_pmatrix, bool covariates, bool get_pval, bool get_connectivity,
             size_t n_voxels, size_t n_subjects, size_t n_permutations, size_t n_covariates, std::vector<int> selected_devices){



  int devices = selected_devices.size();
  if(n_voxels <= BLOCK_SIZE_3){
    batch_size = n_voxels;
    iterations = 1;
    free_voxels = batch_size;

  }else{
    batch_size = BLOCK_SIZE_3;
    iterations = ceil(float(n_voxels)/float(batch_size));
    free_voxels = n_voxels % batch_size;
    if(n_voxels % batch_size == 0)
      free_voxels = batch_size;
  }

  float * h_Y[devices];

  int n_streams = iterations;
  cuda_fphi_variables_per_stream stream_vars[n_streams];
  cuda_fphi_variables_per_device device_vars[devices];

  cudaStream_t device_stream[devices];

  cublasHandle_t device_handles[devices];
  std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

	int blockSize_n_subjects = 32;

	if(n_subjects >= 64){
		blockSize_n_subjects = 64;
	}
	if(n_subjects >= 128){
		blockSize_n_subjects = 128;
	}
	if(n_subjects >= 256){
		blockSize_n_subjects = 256;
	}
	if(n_subjects >= 512){
		blockSize_n_subjects = 512;
	}
	if(n_subjects >= 1024){
		if(n_subjects % 1024 <= n_subjects % 512){
			blockSize_n_subjects = 1024;
		}else{
			blockSize_n_subjects = 512;
		}

	}

	dim3 blockSize(blockSize_n_subjects, 1, 1);
	dim3 gridSize(ceil(float(n_subjects)/float(blockSize_n_subjects)),n_subjects, 1);
/*
	  gpuErrchk(cudaSetDevice(selected_devices[0]));
	  curandState * states;
	  gpuErrchk(cudaMalloc((void**)&states, sizeof(curandState)*n_permutations));
	  unsigned int * pmatrix;
	  gpuErrchk(cudaMalloc((void**)&pmatrix, sizeof(unsigned int)*n_permutations*n_subjects));
	  unsigned int * max;
	  gpuErrchk(cudaMalloc((void**)&max, sizeof(unsigned int)*n_permutations));


		int blockSize_n_permutations = 32;

		if(n_permutations >= 64){
			blockSize_n_permutations = 64;
		}
		if(n_permutations  >= 128){
			blockSize_n_permutations = 128;
		}
		if(n_permutations  >= 256){
			blockSize_n_permutations = 256;
		}
		if(n_permutations  >= 512){
			blockSize_n_permutations = 512;
		}
		if(n_permutations >= 1024){
			blockSize_n_permutations = 1024;
		}

		dim3 blockSize_1(blockSize_n_subjects,1, 1);
		dim3 gridSize_1(ceil(float(n_subjects)/float(blockSize_n_subjects)),n_permutations, 1);


		dim3 blockSize_2(blockSize_n_permutations, 1, 1);
		dim3 gridSize_2(ceil(float(n_permutations)/float(blockSize_n_permutations)), 1, 1);


		init_rand_kernel<<<gridSize_2, blockSize_2, 0, 0>>>(states, time(NULL), n_permutations);
		gpuErrchk(cudaPeekAtLastError());
		init_permutation_matrix<<<gridSize_1,  blockSize_1, 0, 0>>>(pmatrix, max, n_subjects,  n_permutations);
		gpuErrchk(cudaPeekAtLastError());
		create_permutation_matrix<<<gridSize_1,  blockSize_1, 0, 0>>>(states, pmatrix, max,  n_subjects , n_permutations);
		gpuErrchk(cudaPeekAtLastError());
		gpuErrchk(cudaFree(states));
		gpuErrchk(cudaFree(max));

		gpuErrchk(cudaMemcpy(h_pmatrix, pmatrix, sizeof(unsigned int)*n_subjects*n_permutations, cudaMemcpyDeviceToHost));


		for(size_t row = 0 ; row < n_subjects ; row++){
			std::cout << h_pmatrix[row] << "\n";
		}
		std::cin.get();*/

  for(int current_device = 0 ; current_device < devices; current_device++){

	  gpuErrchk(cudaSetDevice(selected_devices[current_device]));
	  device_vars[current_device].device_id = selected_devices[current_device];
	  gpuErrchk(cudaSetDeviceFlags((unsigned int)cudaDeviceScheduleAuto));
	  gpuErrchk(cudaStreamCreate(&device_stream[current_device]));
	  cublasErrchk(cublasCreate_v2(&device_handles[current_device]));
	  cublasErrchk(cublasSetStream_v2(device_handles[current_device], device_stream[current_device]));

  }





  float alpha = 1.f;
  float beta = 0.f;


  float * d_cov[devices];
  float * h_hat = new float[n_subjects*n_subjects];
  if(covariates){
	  compute_hat_CPU(h_cov, h_hat, n_subjects, n_covariates);
  }

  for(int current_device = 0; current_device < devices; current_device++){

	  gpuErrchk(cudaSetDevice(device_vars[current_device].device_id));



	  if(covariates){
		  gpuErrchk(cudaMalloc((void**)&device_vars[current_device].hat, sizeof(float)*n_subjects*n_subjects));
	  }

	  gpuErrchk(cudaMalloc((void**)&device_vars[current_device].evectors, sizeof(float)*n_subjects*n_subjects));



      if(get_pval){
    	  gpuErrchk(cudaMalloc((void**)&device_vars[current_device].pmatrix, sizeof(unsigned int)*n_subjects*n_permutations));
      }

      gpuErrchk(cudaMalloc((void**)&device_vars[current_device].aux_vars.d_Z, sizeof(float)*n_subjects*2));



      gpuErrchk(cudaMalloc((void**)&device_vars[current_device].aux_vars.d_ZTZI, sizeof(float)*2*2));



      if(covariates){


      	   gpuErrchk(cudaMemcpyAsync(device_vars[current_device].hat, h_hat, sizeof(float)*n_subjects*n_subjects,cudaMemcpyHostToDevice, device_stream[current_device]));



      }
      device_vars[current_device].device_id = selected_devices[current_device];
      device_vars[current_device].get_pval = get_pval;
      device_vars[current_device].covariates = covariates;
      device_vars[current_device].get_connectivity = get_connectivity;

      gpuErrchk(cudaMemcpyAsync((float *)device_vars[current_device].aux_vars.d_Z, h_Z,  sizeof(float)*n_subjects*2, cudaMemcpyHostToDevice, device_stream[current_device]));
      if(get_pval)
    	  gpuErrchk(cudaMemcpyAsync((unsigned int *)device_vars[current_device].pmatrix, h_pmatrix, sizeof(unsigned int)*n_subjects*n_permutations, cudaMemcpyHostToDevice, device_stream[current_device]));

      gpuErrchk(cudaMemcpyAsync((float *)device_vars[current_device].evectors, h_evectors, sizeof(float)*n_subjects*n_subjects, cudaMemcpyHostToDevice, device_stream[current_device]));
	  cublasErrchk(cublasSgemm_v2(device_handles[current_device], CUBLAS_OP_T, CUBLAS_OP_N, 2, 2, n_subjects, &alpha, device_vars[current_device].aux_vars.d_Z, n_subjects, device_vars[current_device].aux_vars.d_Z, n_subjects, &beta, device_vars[current_device].aux_vars.d_ZTZI, 2));

	  Inv4by4<<<1, 1, 0, device_stream[current_device]>>>(device_vars[current_device].aux_vars.d_ZTZI);

  }

  for(int current_device = 0 ; current_device < devices; current_device++){

	  gpuErrchk(cudaSetDevice(selected_devices[current_device]));

	  gpuErrchk(cudaStreamSynchronize(device_stream[current_device]));

  }

  stream_vars[0].n_permutations = n_permutations;
  stream_vars[0].n_subjects = n_subjects;
  stream_vars[0].shared_device_vars = device_vars[0];

  stream_vars[0].stream_number= 0;
  stream_vars[0].n_voxels = batch_size;

  size_t total_bytes;
  size_t free_bytes;
  unsigned int total_number_of_threads_usable_by_device[devices];
  unsigned int threads_in_use_per_device[n_streams];
  size_t total_usable_threads = 0;
  size_t used_bytes = run_allocation_test(stream_vars[0]);
  for(int device = 0; device < devices ; device++){
      gpuErrchk(cudaSetDevice(selected_devices[device]));
      gpuErrchk(cudaMemGetInfo(&free_bytes, &total_bytes));

      total_number_of_threads_usable_by_device[device] = floor(float(total_bytes)/float(used_bytes));
      total_usable_threads += total_number_of_threads_usable_by_device[device];
      threads_in_use_per_device[device] = 0;

  }

  int device = 0;
  std::vector<int> current_running_streams;
  pthread_t pthreads[n_streams];
  void * return_val;

  sched_param sched_param;
  sched_param.sched_priority = 3;
  pthread_attr_t attr;
  pthreadErrchk(pthread_attr_init(&attr));
  pthreadErrchk(pthread_attr_setinheritsched(&attr, PTHREAD_EXPLICIT_SCHED));
  pthreadErrchk(pthread_attr_setschedpolicy(&attr, SCHED_RR));
  pthreadErrchk(pthread_attr_setschedparam(&attr, &sched_param));


  if(get_connectivity == true){
	  float * d_sy;
	  float * d_y;
	  float * mean;
	  gpuErrchk(cudaSetDevice(selected_devices[0]));
	  gpuErrchk(cudaMalloc((void**)&d_sy, sizeof(float)*n_subjects*batch_size));
	  gpuErrchk(cudaMalloc((void**)&d_y, sizeof(float)*n_subjects*batch_size));
	  gpuErrchk(cudaMalloc((void**)&mean, sizeof(float)*batch_size));

	  for(size_t stream_number = 0 ; stream_number < iterations ; stream_number++){
		  size_t current_n_voxels;


		  if(stream_number == iterations - 1)
			  current_n_voxels = free_voxels;
		  else
			  current_n_voxels = batch_size;
		  gpuErrchk(cudaMemcpyAsync(d_sy, h_y + stream_number*batch_size*n_subjects, sizeof(float)*n_subjects*current_n_voxels, cudaMemcpyHostToDevice, device_stream[0]));
		  whitened_Y(d_sy, mean,  device_vars[0].evectors,  d_y,  n_subjects,  current_n_voxels, device_handles[0],  device_stream[0]);

		  gpuErrchk(cudaMemcpyAsync(h_y + stream_number*batch_size*n_subjects, d_sy , sizeof(float)*n_subjects*current_n_voxels, cudaMemcpyDeviceToHost, device_stream[0]));

	  }
	  gpuErrchk(cudaStreamSynchronize(device_stream[0]));
	  gpuErrchk(cudaDeviceSynchronize());
	  gpuErrchk(cudaFree(mean));
	  gpuErrchk(cudaFree(d_y));
	  gpuErrchk(cudaFree(d_sy));



	  compute_new_hat_CPU(h_cov, h_hat , h_y,  n_subjects ,  n_covariates,  0);
	  for(int device_index = 0 ; device_index < devices ; device_index++){
		  gpuErrchk(cudaSetDevice(selected_devices[device_index]));
		  gpuErrchk(cudaMemcpy(device_vars[device_index].hat, h_hat, sizeof(float)*n_subjects*n_subjects, cudaMemcpyHostToDevice));
	  }
	  for(size_t current_voxel = 0 ; current_voxel < n_voxels; current_voxel++){


		  compute_new_hat_CPU(h_cov, h_hat , h_y,  n_subjects ,  n_covariates,  current_voxel);
		  for(int device_index = 0 ; device_index < devices ; device_index++){
			  gpuErrchk(cudaSetDevice(selected_devices[device_index]));
			  gpuErrchk(cudaMemcpy(device_vars[device_index].hat, h_hat, sizeof(float)*n_subjects*n_subjects, cudaMemcpyHostToDevice));
		  }

		  for(int current_stream = 0 ; current_stream < n_streams; current_stream++){


			  do{
				  if(threads_in_use_per_device[device] >= total_number_of_threads_usable_by_device[device])
					  device++;
				  else
					  break;
			  }while(device != devices);

			  if(device == devices){
				  int streamIdx;
				  do{

					  for(std::vector<int>::iterator it = current_running_streams.begin(); it != current_running_streams.end() ;it++){
						  pthreadErrchk(pthread_join(pthreads[*it], &return_val));
						  if(return_val != NULL){
							  device = std::distance(selected_devices.begin(), std::find( selected_devices.begin(), selected_devices.end(), stream_vars[*it].shared_device_vars.device_id));
							  streamIdx = *it;
							  threads_in_use_per_device[device] -= 1;
							  break;
						  }
					  }

					  std::this_thread::sleep_for(std::chrono::milliseconds(10));

	  		  }while(device == devices);
				  current_running_streams.erase( std::find(current_running_streams.begin(), current_running_streams.end(), streamIdx));
				  return_val = NULL;
			  }
			  stream_vars[current_stream].shared_device_vars = device_vars[device];

			  stream_vars[current_stream].shared_device_vars.h2 = &conn_h2[current_voxel*n_voxels];
			  if(get_pval)
				  stream_vars[current_stream].shared_device_vars.pvals = &conn_pval[current_voxel*n_voxels];
			  stream_vars[current_stream].shared_device_vars.indicator = &conn_indicator[current_voxel*n_voxels];
			  stream_vars[current_stream].h_y = h_y;
			  stream_vars[current_stream].h_cov = h_cov;

			  stream_vars[current_stream].n_permutations = n_permutations;
			  stream_vars[current_stream].n_subjects = n_subjects;
			  stream_vars[current_stream].n_voxels = batch_size;
			  stream_vars[current_stream].n_covariates = n_covariates;
			  stream_vars[current_stream].stream_number = current_stream;
			  stream_vars[current_stream].cov_number = current_voxel;

			  if(current_stream == iterations - 1)
				  stream_vars[current_stream].n_voxels = free_voxels;
			  pthreadErrchk(pthread_create(&pthreads[current_stream], &attr, run_cudafphi_connect_pthread, (void*)&stream_vars[current_stream]));

			  current_running_streams.push_back(current_stream);
			  threads_in_use_per_device[device] += 1;
			  device++;
			  if(device == devices)
				  device = 0;



		  }

		  for(int current_device = 0 ; current_device < devices; current_device++){
			  threads_in_use_per_device[current_device] = 0;
		  }



		  for(int it = 0; it < current_running_streams.size() ; it++){
			  pthreadErrchk(pthread_join(pthreads[current_running_streams[it]], NULL));
		  }



		  current_running_streams.clear();
		  device = 0;

	  }
  }else{
	  for(int current_stream = 0 ; current_stream < n_streams; current_stream++){


	  	  do{
	  		  if(threads_in_use_per_device[device] >= total_number_of_threads_usable_by_device[device])
	  			  device++;
	  		  else
	  			  break;
	  	  }while(device != devices);

	  	  if(device == devices){
	  		  int streamIdx;
	  		  do{

	  			  for(std::vector<int>::iterator it = current_running_streams.begin(); it != current_running_streams.end() ;it++){
	  				  pthread_tryjoin_np(pthreads[*it], &return_val);
	  				  if(return_val != NULL){
	  					  device = std::distance(selected_devices.begin(), std::find( selected_devices.begin(), selected_devices.end(), stream_vars[*it].shared_device_vars.device_id));
	  					  streamIdx = *it;
	  					  threads_in_use_per_device[device] -= 1;
	  					  break;
	  				  }
	  			  }

	  		  }while(device == devices);

	  		  current_running_streams.erase( std::find(current_running_streams.begin(), current_running_streams.end(), streamIdx));
	  		  return_val = NULL;
	  	  }
	  	  stream_vars[current_stream].shared_device_vars = device_vars[device];

	  	  stream_vars[current_stream].shared_device_vars.h2 = h2;
	  	  stream_vars[current_stream].shared_device_vars.pvals = pvals;
	  	  stream_vars[current_stream].shared_device_vars.indicator = indicator;
	  	  stream_vars[current_stream].h_y = h_y;

	  	  stream_vars[current_stream].n_permutations = n_permutations;
	  	  stream_vars[current_stream].n_subjects = n_subjects;
	  	  stream_vars[current_stream].n_voxels = batch_size;
	  	  stream_vars[current_stream].stream_number = current_stream;

	  	  if(current_stream == iterations - 1)
	  		  stream_vars[current_stream].n_voxels = free_voxels;
	   	  pthreadErrchk(pthread_create(&pthreads[current_stream], NULL, run_cudafphi_pthread, (void*)&stream_vars[current_stream]));

	  	  current_running_streams.push_back(current_stream);
	  	  threads_in_use_per_device[device] += 1;
	  	  device++;
	  	  if(device == devices)
	  		  device = 0;



	    }



	    for(int it = 0; it < current_running_streams.size() ; it++){
	  	  pthreadErrchk(pthread_join(pthreads[current_running_streams[it]], NULL));
	    }

	    current_running_streams.clear();
	    device = 0;

	    for(int current_device = 0 ; current_device < devices; current_device++){
	  	  threads_in_use_per_device[current_device] = 0;
	    }




  }


  for(std::vector<int>::iterator device_it = selected_devices.begin() ; device_it != selected_devices.end(); device_it++){
	  gpuErrchk(cudaSetDevice(*device_it));
	  const size_t  device_Idx = std::distance(selected_devices.begin(), device_it);

	  gpuErrchk(cudaStreamDestroy(device_stream[device_Idx]));
	  gpuErrchk(cudaDeviceSynchronize());
	  gpuErrchk(cudaFree((void*)device_vars[device_Idx].evectors));
      if(get_pval){
    	  gpuErrchk(cudaFree((void*)device_vars[device_Idx].pmatrix));
      }
      gpuErrchk(cudaFree((void*)device_vars[device_Idx].aux_vars.d_Z));
      gpuErrchk(cudaFree((void*)device_vars[device_Idx].aux_vars.d_ZTZI));
    //  gpuErrchk(cudaFreeHost((void*)h_Y[device_Idx]));

	  cublasErrchk(cublasDestroy_v2(device_handles[device_Idx]));

	  if(get_connectivity || covariates){
		  gpuErrchk(cudaFree((void*)device_vars[device_Idx].hat));

	  }

  }

  delete [] h_hat;


  auto elapsed = std::chrono::high_resolution_clock::now() - start;
  long long seconds = std::chrono::duration_cast<std::chrono::seconds>(elapsed).count();
  pthreadErrchk(pthread_attr_destroy(&attr));
  printf ("It took %u seconds.\n", (((int)(seconds))));

  return 1;

}
