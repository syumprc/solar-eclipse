/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
/*
 * cudaFPHI.cu
 * Primary Author: Brian Donohue
 * Date: 12/05/2015
 * Email: bdono09@gmail.edu
 */

/*
 * To whoever wishes to modify and redistribute this file please write your name the date
 * and an explanation of the modification.
 * Modifications:
 *
 */

#include<iostream>

#include<cstdlib>
#include<time.h>
#include<stdio.h>
#include "cudafPHI.cuh"
#include <cuda_runtime.h>


extern size_t iterations;
extern size_t batch_size;
extern size_t free_voxels;

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


/*
static int calculate_batch_size(size_t n_voxels, size_t n_subjects, size_t * batch_size, unsigned int * iterations, size_t * free_voxels){

  *batch_size = n_voxels;
  *iterations = 1;
  *free_voxels = 0;

  size_t freeMem;
  size_t totalMem;
 // cudaMemGetInfo(&freeMem, &totalMem);

  float data_size = 0.f;

  data_size += float(sizeof(float)*n_subjects);
  data_size += float(sizeof(float)*n_subjects);
  data_size += float(sizeof(float));
  data_size += float(sizeof(float));
  data_size += float(sizeof(bool));

  *batch_size = floor((float(float(freeMem)))/data_size);
  *batch_size = n_voxels;
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



  return 1;

}*/


void * run_compute_F_pthread(void * run_compute_F_args){

  const char * file_name = "cudafPHI.cu";

  cudafPHI_variables * cudafPHI_vars = (cudafPHI_variables *)run_compute_F_args;
  cudaSetDevice(cudafPHI_vars->device);
  printf("berfore]\n");
      printError(cudaMemsetAsync(cudafPHI_vars->d_F, 0, sizeof(float)*batch_size*cudafPHI_vars->n_subjects, cudafPHI_vars->stream), file_name, "cudaMemset-d_F", __LINE__);
      printf("after]\n");

  if(cudafPHI_vars->stream_number == (iterations - 1)){
      printError(cudaMemcpyAsync(cudafPHI_vars->d_sy, cudafPHI_vars->h_y + batch_size*cudafPHI_vars->stream_number*cudafPHI_vars->n_subjects, sizeof(float)*cudafPHI_vars->n_subjects*free_voxels, cudaMemcpyHostToDevice, cudafPHI_vars->stream), file_name,
	               "cudaMemcpy-h_sy to d_sy", __LINE__);

      compute_F(cudafPHI_vars->d_hat, cudafPHI_vars->d_sy, cudafPHI_vars->d_evectors,cudafPHI_vars->d_F,
                cudafPHI_vars->compute_F_vars, cudafPHI_vars->covariates,
                cudafPHI_vars->n_subjects,free_voxels, cudafPHI_vars->stream);

  }else{
      printError(cudaMemcpyAsync(cudafPHI_vars->d_sy, cudafPHI_vars->h_y + batch_size*cudafPHI_vars->stream_number*cudafPHI_vars->n_subjects, sizeof(float)*cudafPHI_vars->n_subjects*batch_size, cudaMemcpyHostToDevice, cudafPHI_vars->stream), file_name,
	                  "cudaMemcpy-h_sy to d_sy", __LINE__);

      compute_F(cudafPHI_vars->d_hat, cudafPHI_vars->d_sy, cudafPHI_vars->d_evectors,cudafPHI_vars->d_F,
                cudafPHI_vars->compute_F_vars, cudafPHI_vars->covariates,
                cudafPHI_vars->n_subjects,batch_size, cudafPHI_vars->stream);
  }


  cudaStreamSynchronize(cudafPHI_vars->stream);

}

void * run_compute_h2_pthread(void * run_compute_h2_args){

  cudafPHI_variables * cudafPHI_vars = (cudafPHI_variables *)run_compute_h2_args;
  cudaSetDevice(cudafPHI_vars->device);

  const char * file_name = "cudafPHI.cu";

  if(cudafPHI_vars->stream_number == (iterations - 1)){
      compute_h2(cudafPHI_vars->d_F, cudafPHI_vars->d_h2, cudafPHI_vars->d_indicator, cudafPHI_vars->d_boolean_score,
                 cudafPHI_vars->compute_h2_vars,  cudafPHI_vars->aux_vars,
                 cudafPHI_vars->n_subjects, free_voxels, cudafPHI_vars->stream);

      printError(cudaMemcpyAsync(&cudafPHI_vars->h2[cudafPHI_vars->stream_number*batch_size], cudafPHI_vars->d_h2, sizeof(float)*free_voxels, cudaMemcpyDeviceToHost, cudafPHI_vars->stream), file_name,
	                 "cudaMemcpy-d_h2 to h_h2", __LINE__);

	//  printError(cudaMemcpyAsync(h_boolean_score, d_boolean_score, sizeof(bool)*batch_size, cudaMemcpyDeviceToHost, stream), file_name,
	//                 "cudaMemcpy-h_boolean to d_boolean", __LINE__);
      printError(cudaMemcpyAsync(&cudafPHI_vars->indicator[cudafPHI_vars->stream_number*batch_size], cudafPHI_vars->d_indicator, sizeof(float)*free_voxels, cudaMemcpyDeviceToHost, cudafPHI_vars->stream), file_name,
	                 "cudaMemcpy-d_indicator to h_indicator", __LINE__);

  }else{
      compute_h2(cudafPHI_vars->d_F, cudafPHI_vars->d_h2, cudafPHI_vars->d_indicator, cudafPHI_vars->d_boolean_score,
                 cudafPHI_vars->compute_h2_vars,  cudafPHI_vars->aux_vars,
                 cudafPHI_vars->n_subjects, batch_size, cudafPHI_vars->stream);

      printError(cudaMemcpyAsync(&cudafPHI_vars->h2[cudafPHI_vars->stream_number*batch_size], cudafPHI_vars->d_h2, sizeof(float)*batch_size, cudaMemcpyDeviceToHost, cudafPHI_vars->stream), file_name,
	                 "cudaMemcpy-d_h2 to h_h2", __LINE__);

	//  printError(cudaMemcpyAsync(h_boolean_score, d_boolean_score, sizeof(bool)*batch_size, cudaMemcpyDeviceToHost, stream), file_name,
	//                 "cudaMemcpy-h_boolean to d_boolean", __LINE__);
      printError(cudaMemcpyAsync(&cudafPHI_vars->indicator[cudafPHI_vars->stream_number*batch_size], cudafPHI_vars->d_indicator, sizeof(float)*batch_size, cudaMemcpyDeviceToHost, cudafPHI_vars->stream), file_name,
	                 "cudaMemcpy-d_indicator to h_indicator", __LINE__);
  }

  cudaStreamSynchronize(cudafPHI_vars->stream);

}


void * run_compute_pval_pthread(void * run_compute_pval_args){

  cudafPHI_variables * cudafPHI_vars = (cudafPHI_variables *)run_compute_pval_args;
  cudaSetDevice(cudafPHI_vars->device);

  const char * file_name = "cudafPHI.cu";
  if (cudafPHI_vars->stream_number == (iterations - 1)){


      printError(cudaMemsetAsync(cudafPHI_vars->d_pvals, 0, sizeof(float)*batch_size, cudafPHI_vars->stream),
    		             file_name, "cudaMemset-d_F", __LINE__);
      compute_pvals(cudafPHI_vars->d_sy, cudafPHI_vars->d_F, cudafPHI_vars->d_hat, cudafPHI_vars->compute_h2_vars.d_Sigma_P,
                    cudafPHI_vars->compute_h2_vars.d_Sigma_P, cudafPHI_vars->d_pvals,  cudafPHI_vars->d_pmatrix,
                cudafPHI_vars->d_boolean_score,  cudafPHI_vars->aux_vars,  cudafPHI_vars->pval_vars,  cudafPHI_vars->covariates,
                cudafPHI_vars->n_subjects,  free_voxels, cudafPHI_vars->n_permutations, cudafPHI_vars->stream);
      printError(cudaMemcpyAsync(&cudafPHI_vars->pvals[cudafPHI_vars->stream_number*batch_size], cudafPHI_vars->d_pvals, sizeof(float)*free_voxels, cudaMemcpyDeviceToHost, cudafPHI_vars->stream), file_name,
                 "cudaMemcpy-d_pvals to h_pvals", __LINE__);


  }else{



      printError(cudaMemsetAsync(cudafPHI_vars->d_pvals, 0, sizeof(float)*batch_size, cudafPHI_vars->stream),
    		             file_name, "cudaMemset-d_F", __LINE__);
      compute_pvals(cudafPHI_vars->d_sy, cudafPHI_vars->d_F, cudafPHI_vars->d_hat, cudafPHI_vars->compute_h2_vars.d_Sigma_P,
                    cudafPHI_vars->compute_h2_vars.d_Sigma_P, cudafPHI_vars->d_pvals,  cudafPHI_vars->d_pmatrix,
                cudafPHI_vars->d_boolean_score,  cudafPHI_vars->aux_vars,  cudafPHI_vars->pval_vars,  cudafPHI_vars->covariates,
                cudafPHI_vars->n_subjects,  batch_size, cudafPHI_vars->n_permutations, cudafPHI_vars->stream);
      printError(cudaMemcpyAsync(&cudafPHI_vars->pvals[cudafPHI_vars->stream_number*batch_size], cudafPHI_vars->d_pvals, sizeof(float)*batch_size, cudaMemcpyDeviceToHost, cudafPHI_vars->stream), file_name,
                 "cudaMemcpy-d_pvals to h_pvals", __LINE__);

  }



  cudaStreamSynchronize(cudafPHI_vars->stream);



}


static void cudafPHI(float * h_y, float * h2, float * indicator, float * pvals,float * d_sy,float * d_F,float * d_hat,
             float * d_evectors, float * d_h2, float * d_indicator, float * d_pvals,  unsigned int * d_pmatrix, compute_F_variables compute_F_vars,
             aux_variables aux_vars, compute_h2_variables compute_h2_vars, pval_variables pval_vars,
             bool covariates, bool get_pval, bool * d_boolean_score, size_t n_voxels, size_t n_subjects, size_t n_permutations,
             int current_stream_number, int device, cudaStream_t stream){

 // cudaSetDevice(gpu);
  const char * file_name = "cudafPHI.cu";




 // cudaSetDevice(gpu);

  if(current_stream_number == (iterations - 1)){
      printError(cudaMemcpyAsync(d_sy, h_y + batch_size*current_stream_number*n_subjects, sizeof(float)*n_subjects*free_voxels, cudaMemcpyHostToDevice, stream), file_name,
	               "cudaMemcpy-h_sy to d_sy", __LINE__);

  }else{
      printError(cudaMemcpyAsync(d_sy, h_y + batch_size*current_stream_number*n_subjects, sizeof(float)*n_subjects*batch_size, cudaMemcpyHostToDevice, stream), file_name,
	                  "cudaMemcpy-h_sy to d_sy", __LINE__);
  }

     // cudaDeviceSynchronize();
 //     cudaDeviceSynchronize()
  printError(cudaMemsetAsync(d_F, 0, sizeof(float)*batch_size*n_subjects, stream), file_name, "cudaMemset-d_F", __LINE__);

  if(current_stream_number == (iterations - 1)){
      compute_F(d_hat, d_sy, d_evectors,d_F,
	            compute_F_vars, covariates,
	            n_subjects,free_voxels, stream);

	  compute_h2(d_F, d_h2, d_indicator, d_boolean_score,
	              compute_h2_vars,  aux_vars,
	             n_subjects, free_voxels, stream);


	  printError(cudaMemcpyAsync(&h2[current_stream_number*batch_size], d_h2, sizeof(float)*free_voxels, cudaMemcpyDeviceToHost, stream), file_name,
	                 "cudaMemcpy-d_h2 to h_h2", __LINE__);

	//  printError(cudaMemcpyAsync(h_boolean_score, d_boolean_score, sizeof(bool)*batch_size, cudaMemcpyDeviceToHost, stream), file_name,
	//                 "cudaMemcpy-h_boolean to d_boolean", __LINE__);

	  printError(cudaMemcpyAsync(&indicator[current_stream_number*batch_size], d_indicator, sizeof(float)*free_voxels, cudaMemcpyDeviceToHost, stream), file_name,
	                 "cudaMemcpy-d_indicator to h_indicator", __LINE__);

	//  cudaStreamSynchronize(stream);




	  if(get_pval){

	      printError(cudaMemsetAsync(d_pvals, 0, sizeof(float)*batch_size, stream),
			             file_name, "cudaMemset-d_F", __LINE__);
	      compute_pvals(d_sy, d_F, d_hat, compute_h2_vars.d_Sigma_P,
	                compute_h2_vars.d_Sigma_P, d_pvals,  d_pmatrix,
	                    d_boolean_score,  aux_vars,  pval_vars,  covariates,
	                     n_subjects,  free_voxels, n_permutations, stream);
	      printError(cudaMemcpyAsync(&pvals[current_stream_number*batch_size], d_pvals, sizeof(float)*free_voxels, cudaMemcpyDeviceToHost, stream), file_name,
	                 "cudaMemcpy-d_pvals to h_pvals", __LINE__);
	    //  cudaStreamSynchronize(stream);
	  }




      }else{
	  compute_F(d_hat, d_sy, d_evectors,d_F,
	            compute_F_vars, covariates, n_subjects,
	               batch_size, stream);





	  compute_h2(d_F, d_h2, d_indicator, d_boolean_score,
	              compute_h2_vars,  aux_vars,
	             n_subjects, batch_size, stream);


	  printError(cudaMemcpyAsync(&h2[current_stream_number*batch_size], d_h2, sizeof(float)*batch_size, cudaMemcpyDeviceToHost, stream), file_name,
	                 "cudaMemcpy-d_h2 to h_h2", __LINE__);

	  //printError(cudaMemcpyAsync(h_boolean_score, d_boolean_score, sizeof(bool)*batch_size, cudaMemcpyDeviceToHost, stream), file_name,
	   //             "cudaMemcpy-h_boolean to d_boolean", __LINE__);

	  printError(cudaMemcpyAsync(&indicator[current_stream_number*batch_size], d_indicator, sizeof(float)*batch_size, cudaMemcpyDeviceToHost, stream), file_name,
	                 "cudaMemcpy-d_indicator to h_indicator", __LINE__);

	//  cudaStreamSynchronize(stream);

	  if(get_pval){


	      printError(cudaMemsetAsync(d_pvals, 0, sizeof(float)*batch_size, stream),
			 file_name, "cudaMemset-d_F", __LINE__);

	      compute_pvals(d_sy, d_F, d_hat, compute_h2_vars.d_Sigma_P,
	                compute_h2_vars.d_Sigma_P, d_pvals,  d_pmatrix,
	                    d_boolean_score,  aux_vars,  pval_vars, covariates,
	                     n_subjects,  batch_size, n_permutations, stream);

	      printError(cudaMemcpyAsync(&pvals[current_stream_number*batch_size], d_pvals, sizeof(float)*batch_size, cudaMemcpyDeviceToHost, stream), file_name,
	                 "cudaMemcpy-d_pvals to h_pvals", __LINE__);
	     // cudaStreamSynchronize(stream);
	  }
	 //



      }

    //  cudaStreamSynchronize(stream);











//  printf ("It took %f seconds.\n", (((float)(clock()-timer))/CLOCKS_PER_SEC));
  //for(unsigned int i = 0 ; i < n_voxels ; i++){
//	      std::cout << indicator[i] << " " << h2[i] << " " << pvals[i] << "\n";
	     // h2[iteration*batch_size + i] = h_h2[i];
	     // indicator[iteration*batch_size + i] = h_h2[i];
	     // if(get_pval == true)
		//pvals[iteration*batch_size + i] = h_pvals[i]}
//  }

 // cudaStreamDestroy(stream);

//  cudaDeviceReset();

}

void run_cudafPHI_loop(void * cudafPHI_args){


  const char * file_name = "cudafPHI.cu";
  cudafPHI_variables * cudafPHI_vars;
  cudafPHI_vars = (cudafPHI_variables *)cudafPHI_args;


  printError(cudaMalloc((void**)&cudafPHI_vars->aux_vars.d_Z, sizeof(float)*2*cudafPHI_vars->n_subjects), file_name, "cudaMalloc-d_Z",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->aux_vars.d_ZTZI, sizeof(float)*4), file_name, "cudaMalloc-d_ZTZI",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_F_vars.d_Y, sizeof(float)*cudafPHI_vars->n_subjects), file_name, "cudaMalloc-d_Y",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_F_vars.mean, sizeof(float)), file_name, "cudaMalloc-mean",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_A, sizeof(float)), file_name, "cudaMalloc-d_A",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_B, sizeof(float)), file_name, "cudaMalloc-d_B",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_C, sizeof(float)), file_name, "cudaMalloc-d_C",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_D, sizeof(float)), file_name, "cudaMalloc-d_D",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_E, sizeof(float)), file_name, "cudaMalloc-d_E",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_Sigma_P, sizeof(float)), file_name, "cudaMalloc-d_Sigma_P",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_score, sizeof(float)), file_name, "cudaMalloc-d_score",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_theta, sizeof(float)*2), file_name, "cudaMalloc-d_theta",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_weights, sizeof(float)*cudafPHI_vars->n_subjects), file_name, "cudaMalloc-d_weights",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->pval_vars.syP, sizeof(float)*cudafPHI_vars->n_subjects*(cudafPHI_vars->n_permutations + 1)), file_name, "cudaMalloc-d_syP",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->pval_vars.d_F, sizeof(float)*cudafPHI_vars->n_subjects*(cudafPHI_vars->n_permutations + 1)), file_name, "cudaMalloc-d_F",
         __LINE__);

  printError(cudaMalloc((void**) &cudafPHI_vars->pval_vars.d_score, sizeof(float) * (cudafPHI_vars->n_permutations + 1)),
         file_name, "cudaMalloc-d_score", __LINE__);

  printError(cudaMalloc((void**) &cudafPHI_vars->pval_vars.d_Ts, sizeof(float) * (cudafPHI_vars->n_permutations + 1)),
         file_name, "cudaMalloc-d_Ts", __LINE__);

  printError(cudaMalloc((void**) &cudafPHI_vars->pval_vars.d_a, sizeof(float) * (cudafPHI_vars->n_permutations + 1) * 2),
         file_name, "cudaMalloc-d_a", __LINE__);

  printError(cudaMalloc((void**) &cudafPHI_vars->pval_vars.d_sigmaP, sizeof(float)*(cudafPHI_vars->n_permutations + 1)),
         file_name, "cudaMalloc-d_sigmaP", __LINE__);


  printError(cudaMalloc((void**)&cudafPHI_vars->d_hat, sizeof(float)*cudafPHI_vars->n_subjects*cudafPHI_vars->n_subjects), file_name, "cudaMalloc-d_hat",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->d_evectors, sizeof(float)*cudafPHI_vars->n_subjects*cudafPHI_vars->n_subjects), file_name, "cudaMalloc-d_evectors",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->d_pmatrix, sizeof(float)*cudafPHI_vars->n_subjects*(cudafPHI_vars->n_permutations + 1)), file_name, "cudaMalloc-d_pmatrix",
         __LINE__);



  printError(cudaMalloc((void**)&cudafPHI_vars->d_sy, sizeof(float)*cudafPHI_vars->n_subjects*batch_size),file_name, "cudaMalloc-d_sy",
         __LINE__);


  printError(cudaMalloc((void**)&cudafPHI_vars->d_F, sizeof(float)*batch_size*cudafPHI_vars->n_subjects),
		             file_name, "cudaMalloc-d_F", __LINE__);


  printError(cudaMalloc((void**)&cudafPHI_vars->d_indicator, sizeof(float)*batch_size),
		             file_name, "cudaMalloc-d_indicator", __LINE__);



  printError(cudaMalloc((void**)&cudafPHI_vars->d_h2, sizeof(float)*batch_size),
		             file_name, "cudaMalloc-d_F", __LINE__);



  printError(cudaMalloc((void**)&cudafPHI_vars->d_pvals, sizeof(float)*batch_size),
		             file_name, "cudaMalloc-d_pvals", __LINE__);




  printError(cudaMalloc((void**)&cudafPHI_vars->d_boolean_score, sizeof(bool)*batch_size), file_name, "cudaMalloc-d_boolean_score",
           __LINE__);
  cudaStream_t stream;
  cudaStreamCreate(&stream);


  printError(cudaMemcpyAsync(cudafPHI_vars->d_pmatrix, cudafPHI_vars->h_pmatrix, sizeof(float)*cudafPHI_vars->n_subjects*(cudafPHI_vars->n_permutations + 1), cudaMemcpyHostToDevice, stream), file_name,
         "cudaMemcpy-h_pmatrix to d_pmatrix", __LINE__);

  printError(cudaMemcpyAsync(cudafPHI_vars->d_evectors, cudafPHI_vars->h_evectors, sizeof(float)*cudafPHI_vars->n_subjects*cudafPHI_vars->n_subjects, cudaMemcpyHostToDevice, stream), file_name,
         "cudaMemcpy-evectors to d_evectors", __LINE__);

  printError(cudaMemcpyAsync(cudafPHI_vars->d_hat, cudafPHI_vars->h_hat, sizeof(float)*cudafPHI_vars->n_subjects*cudafPHI_vars->n_subjects, cudaMemcpyHostToDevice,stream), file_name,
             "cudaMemcpy-h_hat to d_hat", __LINE__);

  printError(cudaMemcpyAsync(cudafPHI_vars->aux_vars.d_Z, cudafPHI_vars->h_Z, sizeof(float)*cudafPHI_vars->n_subjects*2, cudaMemcpyHostToDevice, stream), file_name,
             "cudaMemcpy-h_Z to d_Z", __LINE__);


  dim3 blockSize(BLOCK_SIZE_1, 1, 1);
  dim3 gridSize(ceil(float(cudafPHI_vars->n_subjects) / float(BLOCK_SIZE_1)), 1, 1);


  printError(cudaMemsetAsync(cudafPHI_vars->aux_vars.d_ZTZI, 0, sizeof(float)*4, stream),
             file_name, "cudaMemset-d_ZTZI", __LINE__);

  calculate_ZTZ<<<gridSize, blockSize, 0, stream>>>(cudafPHI_vars->aux_vars.d_Z, cudafPHI_vars->aux_vars.d_ZTZI, cudafPHI_vars->n_subjects);

  Inv4by4<<<1, 1, 0, stream>>>(cudafPHI_vars->aux_vars.d_ZTZI);

  for(size_t iteration = 0 ; iteration < iterations ; iteration++){
      cudafPHI(cudafPHI_vars->h_y, cudafPHI_vars->h2, cudafPHI_vars->indicator, cudafPHI_vars->pvals,cudafPHI_vars->d_sy, cudafPHI_vars->d_F,cudafPHI_vars->d_hat,
               cudafPHI_vars->d_evectors, cudafPHI_vars->d_h2, cudafPHI_vars->d_indicator, cudafPHI_vars->d_pvals,  cudafPHI_vars->d_pmatrix, cudafPHI_vars->compute_F_vars,
               cudafPHI_vars->aux_vars, cudafPHI_vars->compute_h2_vars, cudafPHI_vars->pval_vars,
              cudafPHI_vars->covariates, cudafPHI_vars->get_pval, cudafPHI_vars->d_boolean_score,
                 cudafPHI_vars->n_voxels, cudafPHI_vars->n_subjects, cudafPHI_vars->n_permutations,
                iteration, cudafPHI_vars->device, stream);

      cudaStreamSynchronize(stream);
  }


  cudaStreamDestroy(stream);
  cudaDeviceReset();

}

void * run_cudafPHI_thread(void * cudafPHI_args){

}


void * run_allocation_test(void * cudafPHI_args){
  const char * file_name = "cudafPHI.cu";
  cudafPHI_variables * cudafPHI_vars;
  cudafPHI_vars = (cudafPHI_variables *)cudafPHI_args;
  cudaSetDevice(cudafPHI_vars->device);

  printError(cudaMalloc((void**)&cudafPHI_vars->aux_vars.d_Z, sizeof(float)*2*cudafPHI_vars->n_subjects), file_name, "cudaMalloc-d_Z",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->aux_vars.d_ZTZI, sizeof(float)*4), file_name, "cudaMalloc-d_ZTZI",   __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_F_vars.d_Y, sizeof(float)*cudafPHI_vars->n_subjects*batch_size), file_name, "cudaMalloc-d_Y",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_F_vars.mean, sizeof(float)*batch_size), file_name, "cudaMalloc-mean",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_A, sizeof(float)), file_name, "cudaMalloc-d_A",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_B, sizeof(float)), file_name, "cudaMalloc-d_B",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_C, sizeof(float)), file_name, "cudaMalloc-d_C",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_D, sizeof(float)), file_name, "cudaMalloc-d_D",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_E, sizeof(float)), file_name, "cudaMalloc-d_E",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_Sigma_P, sizeof(float)), file_name, "cudaMalloc-d_Sigma_P",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_score, sizeof(float)), file_name, "cudaMalloc-d_score",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_theta, sizeof(float)*2), file_name, "cudaMalloc-d_theta",
         __LINE__);

  printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_weights, sizeof(float)*cudafPHI_vars->n_subjects), file_name, "cudaMalloc-d_weights",
         __LINE__);

  if(cudafPHI_vars->get_pval){
      printError(cudaMalloc((void**)&cudafPHI_vars->pval_vars.syP, sizeof(float)*cudafPHI_vars->n_subjects*(cudafPHI_vars->n_permutations + 1)), file_name, "cudaMalloc-d_syP",
         __LINE__);

      printError(cudaMalloc((void**)&cudafPHI_vars->pval_vars.d_F, sizeof(float)*cudafPHI_vars->n_subjects*(cudafPHI_vars->n_permutations + 1)), file_name, "cudaMalloc-d_F",
         __LINE__);

      printError(cudaMalloc((void**) &cudafPHI_vars->pval_vars.d_score, sizeof(float) * (cudafPHI_vars->n_permutations + 1)),
         file_name, "cudaMalloc-d_score", __LINE__);

      printError(cudaMalloc((void**) &cudafPHI_vars->pval_vars.d_Ts, sizeof(float) * (cudafPHI_vars->n_permutations + 1)),
         file_name, "cudaMalloc-d_Ts", __LINE__);

      printError(cudaMalloc((void**) &cudafPHI_vars->pval_vars.d_a, sizeof(float) * (cudafPHI_vars->n_permutations + 1) * 2),
         file_name, "cudaMalloc-d_a", __LINE__);

      printError(cudaMalloc((void**) &cudafPHI_vars->pval_vars.d_sigmaP, sizeof(float)*(cudafPHI_vars->n_permutations + 1)),
         file_name, "cudaMalloc-d_sigmaP", __LINE__);
      printError(cudaMalloc((void**)&cudafPHI_vars->d_pmatrix, sizeof(float)*cudafPHI_vars->n_subjects*(cudafPHI_vars->n_permutations + 1)), file_name, "cudaMalloc-d_pmatrix",
             __LINE__);
  }

  if(cudafPHI_vars->covariates){
      printError(cudaMalloc((void**)&cudafPHI_vars->d_hat, sizeof(float)*cudafPHI_vars->n_subjects*cudafPHI_vars->n_subjects), file_name, "cudaMalloc-d_hat",
         __LINE__);
  }

  printError(cudaMalloc((void**)&cudafPHI_vars->d_evectors, sizeof(float)*cudafPHI_vars->n_subjects*cudafPHI_vars->n_subjects), file_name, "cudaMalloc-d_evectors",
         __LINE__);






  printError(cudaMalloc((void**)&cudafPHI_vars->d_sy, sizeof(float)*cudafPHI_vars->n_subjects*batch_size),file_name, "cudaMalloc-d_sy",
         __LINE__);


  printError(cudaMalloc((void**)&cudafPHI_vars->d_F, sizeof(float)*batch_size*cudafPHI_vars->n_subjects),
		             file_name, "cudaMalloc-d_F", __LINE__);


  printError(cudaMalloc((void**)&cudafPHI_vars->d_indicator, sizeof(float)*batch_size),
		             file_name, "cudaMalloc-d_indicator", __LINE__);



  printError(cudaMalloc((void**)&cudafPHI_vars->d_h2, sizeof(float)*batch_size),
		             file_name, "cudaMalloc-d_F", __LINE__);


  if(cudafPHI_vars->get_pval){
    printError(cudaMalloc((void**)&cudafPHI_vars->d_pvals, sizeof(float)*batch_size),
		             file_name, "cudaMalloc-d_pvals", __LINE__);
  }




  printError(cudaMalloc((void**)&cudafPHI_vars->d_boolean_score, sizeof(bool)*batch_size), file_name, "cudaMalloc-d_boolean_score",
           __LINE__);


}

void * run_cudafPHI_pthread(void * cudafPHI_args){

    const char * file_name = "cudafPHI.cu";
    cudafPHI_variables * cudafPHI_vars;
    cudafPHI_vars = (cudafPHI_variables *)cudafPHI_args;
    cudaSetDevice(cudafPHI_vars->device);
    printError(cudaMalloc((void**)&cudafPHI_vars->aux_vars.d_Z, sizeof(float)*2*cudafPHI_vars->n_subjects), file_name, "cudaMalloc-d_Z",
           __LINE__);

    printError(cudaMalloc((void**)&cudafPHI_vars->aux_vars.d_ZTZI, sizeof(float)*4), file_name, "cudaMalloc-d_ZTZI",   __LINE__);

    printError(cudaMalloc((void**)&cudafPHI_vars->compute_F_vars.d_Y, sizeof(float)*cudafPHI_vars->n_subjects*batch_size), file_name, "cudaMalloc-d_Y",
           __LINE__);

    printError(cudaMalloc((void**)&cudafPHI_vars->compute_F_vars.mean, sizeof(float)*batch_size), file_name, "cudaMalloc-mean",
           __LINE__);

    printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_A, sizeof(float)), file_name, "cudaMalloc-d_A",
           __LINE__);

    printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_B, sizeof(float)), file_name, "cudaMalloc-d_B",
           __LINE__);

    printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_C, sizeof(float)), file_name, "cudaMalloc-d_C",
           __LINE__);

    printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_D, sizeof(float)), file_name, "cudaMalloc-d_D",
           __LINE__);

    printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_E, sizeof(float)), file_name, "cudaMalloc-d_E",
           __LINE__);

    printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_Sigma_P, sizeof(float)), file_name, "cudaMalloc-d_Sigma_P",
           __LINE__);

    printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_score, sizeof(float)), file_name, "cudaMalloc-d_score",
           __LINE__);

    printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_theta, sizeof(float)*2), file_name, "cudaMalloc-d_theta",
           __LINE__);

    printError(cudaMalloc((void**)&cudafPHI_vars->compute_h2_vars.d_weights, sizeof(float)*cudafPHI_vars->n_subjects), file_name, "cudaMalloc-d_weights",
           __LINE__);


    if(cudafPHI_vars->get_pval == true){
        printError(cudaMalloc((void**)&cudafPHI_vars->pval_vars.syP, sizeof(float)*cudafPHI_vars->n_subjects*(cudafPHI_vars->n_permutations + 1)), file_name, "cudaMalloc-d_syP",
           __LINE__);

        printError(cudaMalloc((void**)&cudafPHI_vars->pval_vars.d_F, sizeof(float)*cudafPHI_vars->n_subjects*(cudafPHI_vars->n_permutations + 1)), file_name, "cudaMalloc-d_F",
           __LINE__);

        printError(cudaMalloc((void**) &cudafPHI_vars->pval_vars.d_score, sizeof(float) * (cudafPHI_vars->n_permutations + 1)),
           file_name, "cudaMalloc-d_score", __LINE__);

        printError(cudaMalloc((void**) &cudafPHI_vars->pval_vars.d_Ts, sizeof(float) * (cudafPHI_vars->n_permutations + 1)),
           file_name, "cudaMalloc-d_Ts", __LINE__);

        printError(cudaMalloc((void**) &cudafPHI_vars->pval_vars.d_a, sizeof(float) * (cudafPHI_vars->n_permutations + 1) * 2),
           file_name, "cudaMalloc-d_a", __LINE__);

        printError(cudaMalloc((void**) &cudafPHI_vars->pval_vars.d_sigmaP, sizeof(float)*(cudafPHI_vars->n_permutations + 1)),
           file_name, "cudaMalloc-d_sigmaP", __LINE__);
        printError(cudaMalloc((void**)&cudafPHI_vars->d_pmatrix, sizeof(float)*cudafPHI_vars->n_subjects*(cudafPHI_vars->n_permutations + 1)), file_name, "cudaMalloc-d_pmatrix",
               __LINE__);
    }

    if(cudafPHI_vars->covariates == true){
	printError(cudaMalloc((void**)&cudafPHI_vars->d_hat, sizeof(float)*cudafPHI_vars->n_subjects*cudafPHI_vars->n_subjects), file_name, "cudaMalloc-d_hat",
           __LINE__);
    }

    printError(cudaMalloc((void**)&cudafPHI_vars->d_evectors, sizeof(float)*cudafPHI_vars->n_subjects*cudafPHI_vars->n_subjects), file_name, "cudaMalloc-d_evectors",
           __LINE__);






    printError(cudaMalloc((void**)&cudafPHI_vars->d_sy, sizeof(float)*cudafPHI_vars->n_subjects*batch_size),file_name, "cudaMalloc-d_sy",
           __LINE__);


    printError(cudaMalloc((void**)&cudafPHI_vars->d_F, sizeof(float)*batch_size*cudafPHI_vars->n_subjects),
		             file_name, "cudaMalloc-d_F", __LINE__);


    printError(cudaMalloc((void**)&cudafPHI_vars->d_indicator, sizeof(float)*batch_size),
		             file_name, "cudaMalloc-d_indicator", __LINE__);



    printError(cudaMalloc((void**)&cudafPHI_vars->d_h2, sizeof(float)*batch_size),
		             file_name, "cudaMalloc-d_F", __LINE__);



    if(cudafPHI_vars->get_pval == true){
      printError(cudaMalloc((void**)&cudafPHI_vars->d_pvals, sizeof(float)*batch_size),
  		             file_name, "cudaMalloc-d_pvals", __LINE__);
    }



    printError(cudaMalloc((void**)&cudafPHI_vars->d_boolean_score, sizeof(bool)*batch_size), file_name, "cudaMalloc-d_boolean_score",
             __LINE__);








    cudaStream_t stream;
    cudaStreamCreate(&cudafPHI_vars->stream);

    if(cudafPHI_vars->get_pval == true){
	printError(cudaMemcpyAsync(cudafPHI_vars->d_pmatrix, cudafPHI_vars->h_pmatrix, sizeof(float)*cudafPHI_vars->n_subjects*(cudafPHI_vars->n_permutations + 1), cudaMemcpyHostToDevice, stream), file_name,
           "cudaMemcpy-h_pmatrix to d_pmatrix", __LINE__);
    }

    printError(cudaMemcpyAsync(cudafPHI_vars->d_evectors, cudafPHI_vars->h_evectors, sizeof(float)*cudafPHI_vars->n_subjects*cudafPHI_vars->n_subjects, cudaMemcpyHostToDevice, stream), file_name,
           "cudaMemcpy-evectors to d_evectors", __LINE__);

    if(cudafPHI_vars->covariates == true){
	printError(cudaMemcpyAsync(cudafPHI_vars->d_hat, cudafPHI_vars->h_hat, sizeof(float)*cudafPHI_vars->n_subjects*cudafPHI_vars->n_subjects, cudaMemcpyHostToDevice,stream), file_name,
               "cudaMemcpy-h_hat to d_hat", __LINE__);
    }

    printError(cudaMemcpyAsync(cudafPHI_vars->aux_vars.d_Z, cudafPHI_vars->h_Z, sizeof(float)*cudafPHI_vars->n_subjects*2, cudaMemcpyHostToDevice, stream), file_name,
               "cudaMemcpy-h_Z to d_Z", __LINE__);


    dim3 blockSize(BLOCK_SIZE_1, 1, 1);
    dim3 gridSize(ceil(float(cudafPHI_vars->n_subjects) / float(BLOCK_SIZE_1)), 1, 1);


    calculate_ZTZ<<<gridSize, blockSize, 0, stream>>>(cudafPHI_vars->aux_vars.d_Z, cudafPHI_vars->aux_vars.d_ZTZI, cudafPHI_vars->n_subjects);

    Inv4by4<<<1, 1, 0, stream>>>(cudafPHI_vars->aux_vars.d_ZTZI);




    cudafPHI(cudafPHI_vars->h_y, cudafPHI_vars->h2, cudafPHI_vars->indicator, cudafPHI_vars->pvals,cudafPHI_vars->d_sy, cudafPHI_vars->d_F,cudafPHI_vars->d_hat,
             cudafPHI_vars->d_evectors, cudafPHI_vars->d_h2, cudafPHI_vars->d_indicator, cudafPHI_vars->d_pvals,  cudafPHI_vars->d_pmatrix, cudafPHI_vars->compute_F_vars,
             cudafPHI_vars->aux_vars, cudafPHI_vars->compute_h2_vars, cudafPHI_vars->pval_vars,
            cudafPHI_vars->covariates, cudafPHI_vars->get_pval, cudafPHI_vars->d_boolean_score,
               cudafPHI_vars->n_voxels, cudafPHI_vars->n_subjects, cudafPHI_vars->n_permutations,
              cudafPHI_vars->stream_number, cudafPHI_vars->device, cudafPHI_vars->stream);


    cudaStreamDestroy(cudafPHI_vars->stream);



}
