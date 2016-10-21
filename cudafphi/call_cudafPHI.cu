#include "cudafPHI.cuh"
#include <iostream>
#include <cuda_runtime.h>
#include<stdio.h>
#include<time.h>
#include<pthread.h>
static inline void printError(cudaError_t error, const char * functionMess, const char * opMess, int Line) {
	if (error == 0){
		return;
	}else{
		fclose(stderr);
		freopen("cudafPHI.err", "w", stderr);
		std::cout << functionMess << " " << opMess << " Line- " << Line << ": "
				<< cudaGetErrorString(error) << "\n";
		std::cerr << functionMess << " " << opMess << " Line- " << Line << ": "
				<< cudaGetErrorString(error) << "\n";
		fclose(stderr);
		cudaDeviceReset();
		exit(0);
	}


}

static __global__ void calculate_ZTZ(float * d_Z, float * ZTZI, size_t n_subjects) {
	size_t rowIdx = threadIdx.x + blockDim.x * blockIdx.x;
	size_t tIdx = threadIdx.x;
	__shared__ float shared_ZTZI_0[BLOCK_SIZE_1];
	__shared__ float shared_ZTZI_2[BLOCK_SIZE_1];
	__shared__ float shared_ZTZI_3[BLOCK_SIZE_1];

	    if (rowIdx < n_subjects) {
	        shared_ZTZI_0[tIdx] = d_Z[rowIdx] * d_Z[rowIdx];
	        shared_ZTZI_2[tIdx] = d_Z[rowIdx + n_subjects] * d_Z[rowIdx];
	        shared_ZTZI_3[tIdx] = d_Z[rowIdx + n_subjects] * d_Z[rowIdx + n_subjects];
	    }else{
	        shared_ZTZI_0[tIdx] = 0.f;
	        shared_ZTZI_2[tIdx] = 0.f;
	        shared_ZTZI_3[tIdx] = 0.f;
	    }
	    __syncthreads();
	    for(unsigned int stride = BLOCK_SIZE_1/2 ; stride > 0 ; stride >>=1){

		if(threadIdx.x < stride && (rowIdx+stride < n_subjects)){
		    shared_ZTZI_0[tIdx] += shared_ZTZI_0[tIdx + stride];
		    shared_ZTZI_2[tIdx] += shared_ZTZI_2[tIdx + stride];
		    shared_ZTZI_3[tIdx] += shared_ZTZI_3[tIdx + stride];
		}

		__syncthreads();
	    }

	    if((tIdx == 0) && (rowIdx < n_subjects)){
		atomicAdd(&ZTZI[0], shared_ZTZI_0[0]);
		atomicAdd(&ZTZI[1], shared_ZTZI_2[0]);
		atomicAdd(&ZTZI[2], shared_ZTZI_2[0]);
		atomicAdd(&ZTZI[3], shared_ZTZI_3[0]);
	    }


}

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


typedef struct{
  aux_variables aux_vars;
  float * d_hat;
  float * h_hat;
  float  * d_evectors;
  unsigned int * d_pmatrix;
  unsigned int * h_pmatrix;
  float * h_Z;
  float * h_evectors;
  size_t n_subjects;
  size_t n_permutations;
  int device_number;
}shared_objects;

static void * allocate_shared_objects_pthread(void * shared_args){

  shared_objects * shared_objs = (shared_objects *)shared_args;

  const char * file_name = "call_cudafPHI.cu";

  cudaSetDevice(shared_objs->device_number);

  printError(cudaMalloc((void**)&shared_objs->aux_vars.d_Z, sizeof(float)*2*shared_objs->n_subjects), file_name, "cudaMalloc-d_Z",
         __LINE__);

  printError(cudaMalloc((void**)&shared_objs->aux_vars.d_ZTZI, sizeof(float)*4), file_name, "cudaMalloc-d_ZTZI",   __LINE__);

  printError(cudaMalloc((void**)&shared_objs->d_hat, sizeof(float)*shared_objs->n_subjects*shared_objs->n_subjects), file_name, "cudaMalloc-d_hat",
         __LINE__);

  printError(cudaMalloc((void**)&shared_objs->d_evectors, sizeof(float)*shared_objs->n_subjects*shared_objs->n_subjects), file_name, "cudaMalloc-d_evectors",
         __LINE__);

  printError(cudaMalloc((void**)&shared_objs->d_pmatrix, sizeof(float)*shared_objs->n_subjects*(shared_objs->n_permutations + 1)), file_name, "cudaMalloc-d_pmatrix",
         __LINE__);

  printError(cudaMemcpyAsync(shared_objs->d_pmatrix, shared_objs->h_pmatrix, sizeof(float)*shared_objs->n_subjects*(shared_objs->n_permutations + 1), cudaMemcpyHostToDevice, 0), file_name,
         "cudaMemcpy-h_pmatrix to d_pmatrix", __LINE__);

  printError(cudaMemcpyAsync(shared_objs->d_evectors, shared_objs->h_evectors, sizeof(float)*shared_objs->n_subjects*shared_objs->n_subjects, cudaMemcpyHostToDevice, 0), file_name,
         "cudaMemcpy-evectors to d_evectors", __LINE__);

  printError(cudaMemcpyAsync(shared_objs->d_hat, shared_objs->h_hat, sizeof(float)*shared_objs->n_subjects*shared_objs->n_subjects, cudaMemcpyHostToDevice,0), file_name,
             "cudaMemcpy-h_hat to d_hat", __LINE__);

  printError(cudaMemcpyAsync(shared_objs->aux_vars.d_Z, shared_objs->h_Z, sizeof(float)*shared_objs->n_subjects*2, cudaMemcpyHostToDevice, 0), file_name,
             "cudaMemcpy-h_Z to d_Z", __LINE__);


  dim3 blockSize(BLOCK_SIZE_1, 1, 1);
  dim3 gridSize(ceil(float(shared_objs->n_subjects) / float(BLOCK_SIZE_1)), 1, 1);


  printError(cudaMemsetAsync(shared_objs->aux_vars.d_ZTZI, 0, sizeof(float)*4, 0),
             file_name, "cudaMemset-d_ZTZI", __LINE__);

  calculate_ZTZ<<<gridSize, blockSize, 0, 0>>>(shared_objs->aux_vars.d_Z, shared_objs->aux_vars.d_ZTZI, shared_objs->n_subjects);

  Inv4by4<<<1, 1, 0, 0>>>(shared_objs->aux_vars.d_ZTZI);

  cudaStreamSynchronize(0);

}


/*
static int calculate_batch_size(size_t n_voxels){

  batch_size = n_voxels;
  iterations = 1;
  free_voxels = 0;

  size_t freeMem;
  size_t totalMem;
  cudaMemGetInfo(&freeMem, &totalMem);

  float data_size = 0.f;

  data_size += float(sizeof(float)*n_subjects);
  data_size += float(sizeof(float)*n_subjects)*10;
  data_size += float(sizeof(float));
  data_size += float(sizeof(float));
  data_size += float(sizeof(bool));

  batch_size = floor(float(float(freeMem)*0.25f)/data_size);
  batch_size = 5;
  if(*batch_size >= n_voxels){
      *batch_size = n_voxels;
      *iterations = 1;
      *free_voxels = n_voxels;
      return 1;
  }


  if(*batch_size == 0)
    return 0;

  *free_voxels = n_voxels%*batch_size;

  *iterations = ceil(float(n_voxels)/float(*batch_size));

  std::cout << freeMem << " " <<  *batch_size << "\n";

  return 1;

}*/


size_t iterations;
size_t batch_size;
size_t free_voxels;


 extern "C" int call_cudafPHI(float * h2, float * indicator, float * pvals,float * h_y,float * h_Z,float * h_hat,
             float * h_evectors, unsigned int * h_pmatrix, bool covariates, bool get_pval,
             size_t n_voxels, size_t n_subjects, size_t n_permutations, std::vector<int> selected_devices){

  const char * file_name = "call_cudafPHI.cu";
  int devices = selected_devices.size();
  if(n_voxels < 1024){
    batch_size = n_voxels;
    iterations = 1;
    free_voxels = batch_size;

  }else{
    batch_size = 1024;
    iterations = ceil(float(n_voxels)/float(batch_size));
    free_voxels = n_voxels % batch_size;
    if(n_voxels % batch_size == 0)
      free_voxels = batch_size;
  }

 //// free_voxels = batch_size;
//  iterations = 1;

  int n_streams = iterations;
  cudafPHI_variables cudafPHI_vars[n_streams];
 /*
  size_t  n_voxel_size[n_streams];

  n_voxel_size[0] = floor(n_voxels/float(n_streams*devices));
  size_t  voxels_calculated[n_streams];
  voxels_calculated[0] = 0;
  cudafPHI_variables cudafPHI_vars[n_streams];
  cudafPHI_vars[0].n_voxels = n_voxel_size[0];
  cudafPHI_vars[0].voxels_used = voxels_calculated[0];
  for(int current_stream = 1; current_stream < n_streams ; current_stream++){

      if(current_stream == (n_streams - 1)){
	  n_voxel_size[current_stream] = n_voxels - (n_streams - 1)*floor(n_voxels/float(n_streams*devices));
      }else{
	  n_voxel_size[current_stream] = floor(n_voxels/float(n_streams*devices));
      }

      voxels_calculated[current_stream] = (current_stream)*floor(n_voxels/float(n_streams*devices));
      cudafPHI_vars[current_stream].n_voxels = n_voxel_size[current_stream];
      cudafPHI_vars[current_stream].voxels_used = voxels_calculated[current_stream];
  }

  aux_variables aux_vars[devices];
  float * d_hat[devices];
  float * d_pmatrix[devices];
  shared_objects shared_objs[devices];
  pthread_t shared_threads[devices];
  for(int device = 0 ; device < devices ; device++){
      shared_objs[device].h_hat = h_hat;
      shared_objs[device].h_pmatrix = h_pmatrix;
      shared_objs[device].h_Z = h_Z;
      shared_objs[device].n_subjects = n_subjects;
      shared_objs[device].h_evectors = h_evectors;
      shared_objs[device].n_permutations = n_permutations;
      shared_objs[device].device_number = device;
     // cudaSetDevice(shared_objs[device].device_number);


     pthread_create(&shared_threads[device], NULL, allocate_shared_objects_pthread, (void*)&shared_objs[device]);

  }


  for(int device = 0 ; device < devices ; device++){
      pthread_join(shared_threads[device], NULL);
  }*/






  int device = 0;



  for(int current_stream = 0; current_stream < n_streams; current_stream++){


      cudaStreamCreate(&cudafPHI_vars[current_stream].stream);

      cudafPHI_vars[current_stream].covariates = covariates;
      cudafPHI_vars[current_stream].indicator = indicator;
      cudafPHI_vars[current_stream].h_y = h_y;

      cudafPHI_vars[current_stream].get_pval = get_pval;

      cudafPHI_vars[current_stream].n_subjects = n_subjects;

      cudafPHI_vars[current_stream].n_permutations = n_permutations;

      if(get_pval){
	  cudafPHI_vars[current_stream].pvals = pvals;
      }


      cudafPHI_vars[current_stream].h2 = h2;


   //   cudafPHI_vars[current_stream].aux_vars = shared_objs[device].aux_vars;
   //   cudafPHI_vars[current_stream].d_evectors = shared_objs[device].d_evectors;
   //   cudafPHI_vars[current_stream].d_hat = shared_objs[device].d_hat;
   //   cudafPHI_vars[current_stream].d_pmatrix = shared_objs[device].d_pmatrix;

      if(covariates){
	  cudafPHI_vars[current_stream].h_hat = h_hat;
      }
      cudafPHI_vars[current_stream].h_evectors = h_evectors;

      if(get_pval){
	  cudafPHI_vars[current_stream].h_pmatrix = h_pmatrix;
      }

      cudafPHI_vars[current_stream].h_Z = h_Z;

      cudafPHI_vars[current_stream].device = selected_devices[device];

      cudafPHI_vars[current_stream].stream_number = current_stream;



      device++;
      if(device == selected_devices.size())
	device = 0;

  }

  clock_t timer = clock();
  pthread_t run_cudafPHI_threads[n_streams];
  pthread_t run_compute_F_threads[n_streams];
  pthread_t run_compute_h2_threads[n_streams];
  pthread_t run_compute_pval_threads[n_streams];

  run_allocation_test((void*)&cudafPHI_vars[0]);

  size_t free_bytes;
  size_t total_bytes;
  size_t total_bytes_to_be_used = free_bytes*n_streams;
  cudaMemGetInfo(&free_bytes, &total_bytes);

  for(int device = 1 ; device < devices ; device++){
      cudaSetDevice(selected_devices[device]);
      size_t total;
      cudaMemGetInfo(&free_bytes, &total);
      total_bytes += total;
  }


  int outer_iterations = ceil(float(total_bytes_to_be_used)/float(total_bytes));
  int streams_in_iteration = ceil(float(n_streams)/float(outer_iterations));
  cudaSetDevice(selected_devices[0]);
  cudaDeviceReset();
  for(int current_outer_iter = 0 ; current_outer_iter < outer_iterations; current_outer_iter++){
      for(int current_stream = 0 ; current_stream < streams_in_iteration; current_stream++){
	  int streamIdx = current_outer_iter*streams_in_iteration + current_stream;
	  if(streamIdx < n_streams){
	      pthread_create(&run_cudafPHI_threads[streamIdx], NULL, run_cudafPHI_pthread, (void*)&cudafPHI_vars[streamIdx]);
	  }

      }

      for(int current_stream = 0 ; current_stream < streams_in_iteration; current_stream++){
	  int streamIdx = current_outer_iter*streams_in_iteration + current_stream;
	  if(streamIdx < n_streams){
	      pthread_join(run_cudafPHI_threads[streamIdx], NULL);
	  }

      }

      for(int device = 0; device < devices; device++){
	  cudaSetDevice(selected_devices[device]);
	  cudaDeviceReset();
      }

  }


  printf ("It took %f seconds.\n", (((float)(clock()-timer))/CLOCKS_PER_SEC));
  /*for(unsigned int i = 0 ; i < n_voxels; i++){
      std::cout << h2[i] << " " << pvals[i] << " " << indicator[i] << "\n";
  }*/

//  timer = clock();
//  cudaSetDevice(0);
//  run_cudafPHI_loop((void*)&cudafPHI_vars[0]);



//  printf ("It took %f seconds.\n", (((float)(clock()-timer))/CLOCKS_PER_SEC));



  return 1;

}
