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
 * cudafphi.cu
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
#include <fstream>
#include<cstdlib>
#include<time.h>
#include<stdio.h>
#include "cudafphi.cuh"
#include <cuda_runtime.h>

extern size_t iterations;
extern size_t batch_size;
extern size_t free_voxels;


static __global__ void calculate_means(const float * d_sy, float * mean, const size_t n_voxels, const size_t n_subjects){

		unsigned int tIdx = threadIdx.x;

		unsigned int rowIdx = blockIdx.x*blockDim.x + threadIdx.x;
		unsigned int vIdx = blockIdx.y*blockDim.y + threadIdx.y;


		extern __shared__ float shared_mean[];

		if(rowIdx < n_subjects && vIdx < n_voxels)
			shared_mean[tIdx] = d_sy[vIdx*n_subjects + rowIdx];
		else
			shared_mean[tIdx] = 0.f;




		__syncthreads();

		if(blockDim.x >= 1024 && tIdx < 512) shared_mean[tIdx] += shared_mean[tIdx + 512];

		__syncthreads();


		if(blockDim.x >= 512 && tIdx < 256) shared_mean[tIdx] += shared_mean[tIdx + 256];

		__syncthreads();



		if(blockDim.x >= 256 && tIdx < 128) shared_mean[tIdx] += shared_mean[tIdx + 128];

		__syncthreads();

		if(blockDim.x >= 128 && tIdx < 64) shared_mean[tIdx] += shared_mean[tIdx + 64];

		__syncthreads();

		if(blockDim.x >= 64 && tIdx < 32) shared_mean[tIdx] += shared_mean[tIdx + 32];

		__syncthreads();

		if(blockDim.x >= 32 && tIdx < 16) shared_mean[tIdx] += shared_mean[tIdx + 16];

		__syncthreads();

		if(blockDim.x >= 16 && tIdx < 8) shared_mean[tIdx] += shared_mean[tIdx + 8];

		__syncthreads();

		if(blockDim.x >= 8 && tIdx < 4) shared_mean[tIdx] += shared_mean[tIdx + 4];

		__syncthreads();

		if(blockDim.x >= 4 && tIdx < 2) shared_mean[tIdx] += shared_mean[tIdx + 2];

		__syncthreads();

		if(blockDim.x >= 2 && tIdx < 1) shared_mean[tIdx] += shared_mean[tIdx + 1];

		__syncthreads();

		if(tIdx == 0){
			atomicAdd(&mean[vIdx], shared_mean[0]);
		}





}

static __global__ void demean_columns(float * d_Y, const float * d_sy, float * mean, size_t n_voxels, size_t n_subjects){

  const size_t rowIdx = threadIdx.x +blockDim.x*blockIdx.x;
  const size_t voxel = threadIdx.y + blockDim.y*blockIdx.y;


  if(rowIdx < n_subjects && voxel < n_voxels){
	  const float value  =  d_sy[voxel*n_subjects + rowIdx] - mean[voxel]/float(n_subjects);
      d_Y[rowIdx + voxel*n_subjects] = value;
  }
}



static __global__ void calculate_sigma(const float * d_sy, float * sigma,  size_t n_voxels,size_t n_subjects){

	const unsigned int tIdx = threadIdx.x;

	const unsigned int rowIdx = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int vIdx = blockIdx.y*blockDim.y + threadIdx.y;


	extern __shared__ float shared_sigma[];

	if(rowIdx < n_subjects && vIdx < n_voxels){
		float const value = d_sy[vIdx*n_subjects + rowIdx];
		shared_sigma[tIdx] = value*value;
	}else{
		shared_sigma[tIdx] = 0.f;
	}




	__syncthreads();

	if(blockDim.x >= 1024 && tIdx < 512) shared_sigma[tIdx] += shared_sigma[tIdx + 512];

	__syncthreads();


	if(blockDim.x >= 512 && tIdx < 256) shared_sigma[tIdx] += shared_sigma[tIdx + 256];

	__syncthreads();



	if(blockDim.x >= 256 && tIdx < 128) shared_sigma[tIdx] += shared_sigma[tIdx + 128];

	__syncthreads();

	if(blockDim.x >= 128 && tIdx < 64) shared_sigma[tIdx] += shared_sigma[tIdx + 64];

	__syncthreads();

	if(blockDim.x >= 64 && tIdx < 32) shared_sigma[tIdx] += shared_sigma[tIdx + 32];

	__syncthreads();

	if(blockDim.x >= 32 && tIdx < 16) shared_sigma[tIdx] += shared_sigma[tIdx + 16];

	__syncthreads();

	if(blockDim.x >= 16 && tIdx < 8) shared_sigma[tIdx] += shared_sigma[tIdx + 8];

	__syncthreads();

	if(blockDim.x >= 8 && tIdx < 4) shared_sigma[tIdx] += shared_sigma[tIdx + 4];

	__syncthreads();

	if(blockDim.x >= 4 && tIdx < 2) shared_sigma[tIdx] += shared_sigma[tIdx + 2];

	__syncthreads();

	if(blockDim.x >= 2 && tIdx < 1) shared_sigma[tIdx] += shared_sigma[tIdx + 1];

	__syncthreads();


	if(tIdx == 0){
		atomicAdd(&sigma[vIdx], shared_sigma[0]);

	}



}





static __global__ void calculate_inverse_normal(float * d_Y,  const float * d_sigma, size_t n_voxels, size_t n_subjects){

  size_t rowIdx = threadIdx.x +blockDim.x*blockIdx.x;
  size_t voxel = threadIdx.y + blockDim.y*blockIdx.y;


  if(rowIdx < n_subjects && voxel < n_voxels){
	  float value = d_Y[rowIdx + voxel*n_subjects]/sqrt(d_sigma[voxel]);
      d_Y[rowIdx + voxel*n_subjects] = value;
  }
}

void whitened_Y(float * d_sy, float * mean, const float * d_evectors, float * d_y, size_t n_subjects, size_t n_voxels,cublasHandle_t handle, cudaStream_t stream){

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

	dim3 blockSize_mean(blockSize_n_subjects, 1, 1);
	dim3 gridSize_mean(ceil(float(n_subjects)/float(blockSize_n_subjects)), n_voxels, 1);

	dim3 blockSize_set(blockSize_n_subjects, 1024/blockSize_n_subjects, 1);

	dim3 gridSize_set(ceil(float(n_subjects)/float(blockSize_n_subjects)), ceil(float(n_voxels)/float(1024.f/blockSize_n_subjects)), 1);

    gpuErrchk(cudaMemsetAsync(mean, 0, sizeof(float)*n_voxels, stream ));

    calculate_means<<<gridSize_mean, blockSize_mean, sizeof(float)*blockSize_n_subjects, stream>>>(d_sy, mean, n_voxels, n_subjects);
    gpuErrchk(cudaPeekAtLastError());
	demean_columns<<<gridSize_set, blockSize_set, 0, stream>>>(d_y, d_sy, mean, n_voxels, n_subjects);
    gpuErrchk(cudaPeekAtLastError());

    gpuErrchk(cudaMemsetAsync(mean, 0, sizeof(float)*n_voxels, stream ));

    calculate_sigma<<<gridSize_mean, blockSize_mean, sizeof(float)*blockSize_n_subjects, stream>>>(d_y, mean, n_voxels, n_subjects);
    gpuErrchk(cudaPeekAtLastError());


    calculate_inverse_normal<<<gridSize_set, blockSize_set, 0 , stream>>>(d_y, mean, n_voxels, n_subjects);
    gpuErrchk(cudaPeekAtLastError());

	float alpha = 1.f;
	float beta = 0.f;

	cublasErrchk(cublasSgemm_v2(handle, CUBLAS_OP_T, CUBLAS_OP_N, n_subjects,n_voxels, n_subjects, &alpha, d_evectors,
			n_subjects, d_y, n_subjects, &beta, d_sy, n_subjects));



}





static  void cudafphi(float * h2, float * indicator, float * pvals,float * d_sy,float * d_F, const float * d_hat,
              const float * d_evectors, float * d_h2, float * d_indicator, float * d_pvals,    compute_F_variables compute_F_vars,
             aux_variables aux_vars, compute_h2_variables compute_h2_vars, pval_variables pval_vars, const unsigned int * d_pmatrix,
             bool covariates, bool get_pval, bool * d_boolean_score, size_t n_voxels, size_t n_subjects, size_t n_permutations,
             int current_stream_number, cudaStream_t stream, cublasHandle_t handle){

  gpuErrchk(cudaMemsetAsync(d_F, 0, sizeof(float)*n_voxels*n_subjects, stream));
  compute_F(d_hat, d_sy, d_evectors,d_F,
            compute_F_vars, covariates,
            n_subjects,n_voxels, stream, handle);

  compute_h2(d_F, d_h2, d_indicator, d_boolean_score,
          compute_h2_vars,  aux_vars,
             n_subjects, n_voxels, stream);



  gpuErrchk(cudaMemcpyAsync(&h2[current_stream_number*batch_size], d_h2, sizeof(float)*n_voxels, cudaMemcpyDeviceToHost, stream));


  gpuErrchk(cudaMemcpyAsync(&indicator[current_stream_number*batch_size], d_indicator, sizeof(float)*n_voxels, cudaMemcpyDeviceToHost, stream));

  if(get_pval){
      gpuErrchk(cudaMemsetAsync(d_pvals, 0, sizeof(float)*n_voxels, stream));
      compute_pvals(d_sy, d_F, d_hat, compute_h2_vars.d_Sigma_A, compute_h2_vars.d_Sigma_E,
               d_pvals, (const unsigned int *)d_pmatrix,
                    d_boolean_score,  aux_vars,  pval_vars, covariates,
                     n_subjects,  n_voxels, n_permutations, stream, handle);
      gpuErrchk(cudaMemcpyAsync(&pvals[current_stream_number*batch_size], d_pvals, sizeof(float)*n_voxels, cudaMemcpyDeviceToHost, stream));
  }


}


static  void cudafphi_connectivity(float * h2, float * indicator, float * pvals,float * d_sy,float * d_F, const float * d_hat,
              const float * d_evectors, float * d_h2, float * d_indicator, float * d_pvals,    compute_F_variables compute_F_vars,
             aux_variables aux_vars, compute_h2_variables compute_h2_vars, pval_variables pval_vars, const unsigned int * d_pmatrix,
             bool covariates, bool get_pval, bool * d_boolean_score, size_t n_voxels, size_t n_subjects, size_t n_permutations,
             int current_stream_number, cudaStream_t stream, cublasHandle_t handle){



  float alpha = 1.f;
  float beta = 0.f;

  cublasErrchk(cublasSgemm_v2(handle, CUBLAS_OP_N, CUBLAS_OP_N, n_subjects,n_voxels, n_subjects, &alpha, d_hat, n_subjects, d_sy, n_subjects, &beta, d_F, n_subjects));


  compute_h2(d_F, d_h2, d_indicator, d_boolean_score,
          compute_h2_vars,  aux_vars,
             n_subjects, n_voxels, stream);
  gpuErrchk(cudaMemcpyAsync(&h2[current_stream_number*batch_size], d_h2, sizeof(float)*n_voxels, cudaMemcpyDeviceToHost, stream));


  gpuErrchk(cudaMemcpyAsync(&indicator[current_stream_number*batch_size], d_indicator, sizeof(float)*n_voxels, cudaMemcpyDeviceToHost, stream));

  if(get_pval){
      gpuErrchk(cudaMemsetAsync(d_pvals, 0, sizeof(float)*n_voxels, stream));
      compute_pvals(d_sy, d_F, d_hat, compute_h2_vars.d_Sigma_A, compute_h2_vars.d_Sigma_E,
               d_pvals, (const unsigned int *)d_pmatrix,
                    d_boolean_score,  aux_vars,  pval_vars, covariates,
                     n_subjects,  n_voxels, n_permutations, stream, handle);
      gpuErrchk(cudaMemcpyAsync(&pvals[current_stream_number*batch_size], d_pvals, sizeof(float)*n_voxels, cudaMemcpyDeviceToHost, stream));
  }
}


size_t run_allocation_test(cuda_fphi_variables_per_stream cudafphi_vars){


	gpuErrchk(cudaSetDevice(cudafphi_vars.shared_device_vars.device_id));

	compute_h2_variables compute_h2_vars;
	compute_F_variables compute_F_vars;
	pval_variables pval_vars;
	float * d_sy;
	float * d_F;
	bool * d_boolean_score;
	float * d_indicator;
	float * d_h2;
	float * d_pval;
	size_t n_subjects = cudafphi_vars.n_subjects;
	size_t n_permutations = cudafphi_vars.n_permutations;


	gpuErrchk(cudaMalloc((void**)&compute_F_vars.d_Y, sizeof(float)*n_subjects*batch_size));

	gpuErrchk(cudaMalloc((void**)&compute_F_vars.mean_or_sigma, sizeof(float)*batch_size));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_A, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_B, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_C, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_D, sizeof(float)) );

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_E, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_Sigma_P, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_score, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_theta, sizeof(float)*2));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_weights, sizeof(float)*n_subjects));
	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_Sigma_A, sizeof(float)*batch_size));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_Sigma_E, sizeof(float)*batch_size));

	if(cudafphi_vars.shared_device_vars.get_pval){
		gpuErrchk(cudaMalloc((void**)&pval_vars.syP, sizeof(float)*n_subjects*(n_permutations + 1)));

		gpuErrchk(cudaMalloc((void**)&pval_vars.d_F, sizeof(float)*n_subjects*(n_permutations + 1)));

		gpuErrchk(cudaMalloc((void**) &pval_vars.d_score, sizeof(float) * (n_permutations + 1)));

		gpuErrchk(cudaMalloc((void**) &pval_vars.d_Ts, sizeof(float) * (n_permutations + 1)));

		gpuErrchk(cudaMalloc((void**) &pval_vars.d_a, sizeof(float) * (n_permutations + 1) * 2));

    	 gpuErrchk(cudaMalloc((void**) &pval_vars.d_sigmaP, sizeof(float)*(n_permutations + 1)));

	}


	gpuErrchk(cudaMalloc((void**)&d_sy, sizeof(float)*n_subjects*batch_size));


	gpuErrchk(cudaMalloc((void**)&d_F, sizeof(float)*batch_size*n_subjects));


	gpuErrchk(cudaMalloc((void**)&d_indicator, sizeof(float)*batch_size));



	gpuErrchk(cudaMalloc((void**)&d_h2, sizeof(float)*batch_size));


	if(cudafphi_vars.shared_device_vars.get_pval){
		gpuErrchk(cudaMalloc((void**)&d_pval, sizeof(float)*batch_size));
	}




	gpuErrchk(cudaMalloc((void**)&d_boolean_score, sizeof(bool)*batch_size));
	size_t freeMem;
	size_t totalMem;
	size_t usedMem;
	gpuErrchk(cudaMemGetInfo(&freeMem, &totalMem));
	usedMem = totalMem - freeMem;

	gpuErrchk(cudaFree(compute_F_vars.d_Y));


	gpuErrchk(cudaFree(compute_F_vars.mean_or_sigma));



	gpuErrchk(cudaFree(compute_h2_vars.d_A));

	gpuErrchk(cudaFree(compute_h2_vars.d_B));

	gpuErrchk(cudaFree(compute_h2_vars.d_C));

	gpuErrchk(cudaFree(compute_h2_vars.d_D));

	gpuErrchk(cudaFree(compute_h2_vars.d_E));

	gpuErrchk(cudaFree(compute_h2_vars.d_Sigma_P));

	gpuErrchk(cudaFree(compute_h2_vars.d_score));

	gpuErrchk(cudaFree(compute_h2_vars.d_theta));

	gpuErrchk(cudaFree(compute_h2_vars.d_weights));


	gpuErrchk(cudaFree(compute_h2_vars.d_Sigma_E));

	gpuErrchk(cudaFree(compute_h2_vars.d_Sigma_A));

	if(cudafphi_vars.shared_device_vars.get_pval){
		gpuErrchk(cudaFree(pval_vars.syP));

		gpuErrchk(cudaFree(pval_vars.d_F));

		gpuErrchk(cudaFree( pval_vars.d_score));

		gpuErrchk(cudaFree( pval_vars.d_Ts));

		gpuErrchk(cudaFree( pval_vars.d_a));

		gpuErrchk(cudaFree( pval_vars.d_sigmaP));

	}


	gpuErrchk(cudaFree(d_sy));


	gpuErrchk(cudaFree(d_F));


	gpuErrchk(cudaFree(d_indicator));



	gpuErrchk(cudaFree(d_h2));


	gpuErrchk(cudaFree(d_boolean_score));


	if(cudafphi_vars.shared_device_vars.get_pval){
		gpuErrchk(cudaFree(d_pval));
	}

	return  usedMem;

}

void * run_cudafphi_pthread(void * cudafphi_args){

   	cuda_fphi_variables_per_stream * cudafphi_vars;
    cudafphi_vars = (cuda_fphi_variables_per_stream *)cudafphi_args;
	gpuErrchk(cudaSetDevice(cudafphi_vars->shared_device_vars.device_id));

    compute_F_variables compute_F_vars;
    compute_h2_variables compute_h2_vars;
    pval_variables pval_vars;

    float * d_sy;
    float * d_F;
    bool * d_boolean_score;
    float * d_pvals;
    float * d_h2;
    float * d_indicator;


	size_t n_subjects = cudafphi_vars->n_subjects;
	size_t n_permutations = cudafphi_vars->n_permutations;
	size_t n_voxels = cudafphi_vars->n_voxels;



	gpuErrchk(cudaMalloc((void**)&compute_F_vars.d_Y, sizeof(float)*n_subjects*n_voxels));

	gpuErrchk(cudaMalloc((void**)&compute_F_vars.mean_or_sigma, sizeof(float)*n_voxels));
	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_A, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_B, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_C, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_D, sizeof(float)) );

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_E, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_Sigma_P, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_score, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_theta, sizeof(float)*2));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_weights, sizeof(float)*n_subjects));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_Sigma_A, sizeof(float)*n_voxels));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_Sigma_E, sizeof(float)*n_voxels));

	if(cudafphi_vars->shared_device_vars.get_pval){
		gpuErrchk(cudaMalloc((void**)&pval_vars.syP, sizeof(float)*n_subjects*(n_permutations + 1)));

		gpuErrchk(cudaMalloc((void**)&pval_vars.d_F, sizeof(float)*n_subjects*(n_permutations + 1)));

		gpuErrchk(cudaMalloc((void**) &pval_vars.d_score, sizeof(float) * (n_permutations + 1)));

		gpuErrchk(cudaMalloc((void**) &pval_vars.d_Ts, sizeof(float) * (n_permutations + 1)));

		gpuErrchk(cudaMalloc((void**) &pval_vars.d_a, sizeof(float) * (n_permutations + 1) * 2));

    	 gpuErrchk(cudaMalloc((void**) &pval_vars.d_sigmaP, sizeof(float)*(n_permutations + 1)));

	}


	gpuErrchk(cudaMalloc((void**)&d_sy, sizeof(float)*n_subjects*n_voxels));


	gpuErrchk(cudaMalloc((void**)&d_F, sizeof(float)*n_voxels*n_subjects));


	gpuErrchk(cudaMalloc((void**)&d_indicator, sizeof(float)*n_voxels));



	gpuErrchk(cudaMalloc((void**)&d_h2, sizeof(float)*n_voxels));


	if(cudafphi_vars->shared_device_vars.get_pval){
		gpuErrchk(cudaMalloc((void**)&d_pvals, sizeof(float)*n_voxels));
	}




	gpuErrchk(cudaMalloc((void**)&d_boolean_score, sizeof(bool)*n_voxels));



    cudaStream_t stream;

    cublasHandle_t handle;

    cublasErrchk(cublasCreate_v2(&handle));
    gpuErrchk(cudaStreamCreate(&stream));
    cublasErrchk(cublasSetStream_v2(handle, stream));
    float alpha = 1;
    float beta = 0.f;

    gpuErrchk(cudaMemcpyAsync(d_sy, cudafphi_vars->h_y + batch_size*cudafphi_vars->stream_number*cudafphi_vars->n_subjects, sizeof(float)*cudafphi_vars->n_subjects*n_voxels, cudaMemcpyHostToDevice, stream));


    cudafphi(cudafphi_vars->shared_device_vars.h2, cudafphi_vars->shared_device_vars.indicator, cudafphi_vars->shared_device_vars.pvals, d_sy, d_F,cudafphi_vars->shared_device_vars.hat,
    			cudafphi_vars->shared_device_vars.evectors, d_h2, d_indicator, d_pvals,  compute_F_vars,
    			cudafphi_vars->shared_device_vars.aux_vars, compute_h2_vars, pval_vars, cudafphi_vars->shared_device_vars.pmatrix,
    			cudafphi_vars->shared_device_vars.covariates,cudafphi_vars->shared_device_vars.get_pval, d_boolean_score,
    			n_voxels, cudafphi_vars->n_subjects, cudafphi_vars->n_permutations,
    			cudafphi_vars->stream_number, stream, handle);



  //  if(cudafphi_vars->shared_device_vars.get_connectivity){
        gpuErrchk(cudaMemcpyAsync(cudafphi_vars->h_y+ batch_size*cudafphi_vars->stream_number*cudafphi_vars->n_subjects, d_sy, sizeof(float)*cudafphi_vars->n_subjects*n_voxels, cudaMemcpyDeviceToHost, stream));

   // }




     gpuErrchk(cudaStreamSynchronize(stream));


     std::ofstream y_out("Y.mat.csv");

     for(size_t row = 0 ; row < n_subjects; row++){
    	 for(size_t voxel = 0 ; voxel < n_voxels ; voxel++){
    		 if(voxel + 1 == n_voxels)
    			 y_out << cudafphi_vars->h_y[voxel*n_subjects + row] << "\n";
    		 else
    			 y_out << cudafphi_vars->h_y[voxel*n_subjects + row] << ",";
    	 }

     }

     y_out.close();

    gpuErrchk(cudaStreamDestroy(stream));
    cublasErrchk(cublasDestroy_v2(handle));


	gpuErrchk(cudaFree(compute_F_vars.d_Y));


	gpuErrchk(cudaFree(compute_F_vars.mean_or_sigma));



	gpuErrchk(cudaFree(compute_h2_vars.d_A));

	gpuErrchk(cudaFree(compute_h2_vars.d_B));

	gpuErrchk(cudaFree(compute_h2_vars.d_C));

	gpuErrchk(cudaFree(compute_h2_vars.d_D));

	gpuErrchk(cudaFree(compute_h2_vars.d_E));

	gpuErrchk(cudaFree(compute_h2_vars.d_Sigma_P));

	gpuErrchk(cudaFree(compute_h2_vars.d_score));

	gpuErrchk(cudaFree(compute_h2_vars.d_theta));

	gpuErrchk(cudaFree(compute_h2_vars.d_weights));


	gpuErrchk(cudaFree(compute_h2_vars.d_Sigma_E));

	gpuErrchk(cudaFree(compute_h2_vars.d_Sigma_A));

	if(cudafphi_vars->shared_device_vars.get_pval){
		gpuErrchk(cudaFree(pval_vars.syP));

		gpuErrchk(cudaFree(pval_vars.d_F));

		gpuErrchk(cudaFree( pval_vars.d_score));

		gpuErrchk(cudaFree( pval_vars.d_Ts));

		gpuErrchk(cudaFree( pval_vars.d_a));

		gpuErrchk(cudaFree( pval_vars.d_sigmaP));

	}


	gpuErrchk(cudaFree(d_sy));


	gpuErrchk(cudaFree(d_F));


	gpuErrchk(cudaFree(d_indicator));



	gpuErrchk(cudaFree(d_h2));


	gpuErrchk(cudaFree(d_boolean_score));


	if(cudafphi_vars->shared_device_vars.get_pval){
		gpuErrchk(cudaFree(d_pvals));
	}

    int * success = new int;
    *success = 1;
    return (void*)success;

}

static __global__ void  compute_Inverse(float * d_A, size_t n_covariates){
	int i, j, k;
	float ratio;
	float a;
	int n = n_covariates;
    for(i = 0; i < n; i++){
        for(j = n; j < 2*n; j++){
            if(i==(j-n))
                d_A[j*n_covariates +i] = 1.0;
            else
            	d_A[j*n_covariates +i]= 0.0;
        }
    }
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            if(i!=j){
                ratio = d_A[i*n_covariates +j]/ d_A[i*n_covariates +i];
                for(k = 0; k < 2*n; k++){
                	 d_A[k*n_covariates +j] -= ratio *  d_A[k*n_covariates +i];
                }
            }
        }
    }
    for(i = 0; i < n; i++){
        a = d_A[i*n_covariates +i];
        for(j = 0; j < 2*n; j++){
        	d_A[j*n_covariates +i] /= a;
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

static void compute_hat_matrix(const float * d_cov, float * inverse_d_cov, float * d_hat, float * intermediate_hat,  cublasHandle_t handle, cudaStream_t stream, size_t n_subjects, size_t n_covariates){

	float alpha = 1.f;
	float beta = 0.f;

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


	dim3 blockSize_set(blockSize_n_subjects,1 , 1);

	dim3 gridSize_set(ceil(float(n_subjects)/float(blockSize_n_subjects)), n_subjects, 1);



	cublasErrchk(cublasSgemm_v2(handle, CUBLAS_OP_T, CUBLAS_OP_N, n_covariates, n_covariates, n_subjects, &alpha, d_cov, n_subjects, d_cov, n_subjects, &beta, inverse_d_cov, n_covariates));
//	cublasErrchk(cublasSetPointerMode_v2(handle,CUBLAS_POINTER_MODE_HOST));
	compute_Inverse<<<1, 1, 0 , stream>>>(inverse_d_cov, n_covariates);


//	cublasErrchk(cublasSmatinvBatched(handle, n_covariates, (const float**)&inverse_d_cov, n_covariates, (float **)&inverse_d_cov_output, n_covariates, info, 1));

//	cublasErrchk(cublasSetPointerMode_v2(handle,CUBLAS_POINTER_MODE_DEVICE));


	cublasErrchk(cublasSgemm_v2(handle, CUBLAS_OP_N, CUBLAS_OP_N, n_subjects, n_covariates, n_covariates, &alpha, d_cov, n_subjects, inverse_d_cov, n_covariates, &beta, intermediate_hat, n_subjects));

	cublasErrchk(cublasSgemm_v2(handle, CUBLAS_OP_N, CUBLAS_OP_T, n_subjects, n_subjects, n_covariates, &alpha, intermediate_hat, n_subjects, d_cov, n_subjects, &beta, d_hat, n_subjects));

	subtract_matrices<<<gridSize_set, blockSize_set, 0, stream>>>(d_hat, n_subjects);



}


void * run_cudafphi_connect_pthread(void * cudafphi_args){

   	cuda_fphi_variables_per_stream * cudafphi_vars;
    cudafphi_vars = (cuda_fphi_variables_per_stream *)cudafphi_args;
	gpuErrchk(cudaSetDevice(cudafphi_vars->shared_device_vars.device_id));

    compute_F_variables compute_F_vars;
    compute_h2_variables compute_h2_vars;
    pval_variables pval_vars;

    float * d_sy;
    float * d_F;
    bool * d_boolean_score;
    float * d_pvals;
    float * d_h2;
    float * d_indicator;


	size_t n_subjects = cudafphi_vars->n_subjects;
	size_t n_permutations = cudafphi_vars->n_permutations;
	size_t n_voxels = cudafphi_vars->n_voxels;




	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_A, sizeof(float)));



	gpuErrchk(cudaMalloc((void**)&compute_F_vars.d_Y, sizeof(float)*n_subjects*n_voxels));

	gpuErrchk(cudaMalloc((void**)&compute_F_vars.mean_or_sigma, sizeof(float)*n_voxels));
	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_A, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_B, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_C, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_D, sizeof(float)) );

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_E, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_Sigma_P, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_score, sizeof(float)));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_theta, sizeof(float)*2));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_weights, sizeof(float)*n_subjects));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_Sigma_A, sizeof(float)*n_voxels));

	gpuErrchk(cudaMalloc((void**)&compute_h2_vars.d_Sigma_E, sizeof(float)*n_voxels));

	if(cudafphi_vars->shared_device_vars.get_pval){
		gpuErrchk(cudaMalloc((void**)&pval_vars.syP, sizeof(float)*n_subjects*(n_permutations + 1)));

		gpuErrchk(cudaMalloc((void**)&pval_vars.d_F, sizeof(float)*n_subjects*(n_permutations + 1)));

		gpuErrchk(cudaMalloc((void**) &pval_vars.d_score, sizeof(float) * (n_permutations + 1)));

		gpuErrchk(cudaMalloc((void**) &pval_vars.d_Ts, sizeof(float) * (n_permutations + 1)));

		gpuErrchk(cudaMalloc((void**) &pval_vars.d_a, sizeof(float) * (n_permutations + 1) * 2));

    	 gpuErrchk(cudaMalloc((void**) &pval_vars.d_sigmaP, sizeof(float)*(n_permutations + 1)));

	}


	gpuErrchk(cudaMalloc((void**)&d_sy, sizeof(float)*n_subjects*n_voxels));


	gpuErrchk(cudaMalloc((void**)&d_F, sizeof(float)*n_voxels*n_subjects));


	gpuErrchk(cudaMalloc((void**)&d_indicator, sizeof(float)*n_voxels));



	gpuErrchk(cudaMalloc((void**)&d_h2, sizeof(float)*n_voxels));


	if(cudafphi_vars->shared_device_vars.get_pval){
		gpuErrchk(cudaMalloc((void**)&d_pvals, sizeof(float)*n_voxels));
	}




	gpuErrchk(cudaMalloc((void**)&d_boolean_score, sizeof(bool)*n_voxels));



    cudaStream_t stream;

    cublasHandle_t handle;

    cublasErrchk(cublasCreate_v2(&handle));
    gpuErrchk(cudaStreamCreate(&stream));
    cublasErrchk(cublasSetStream_v2(handle, stream));
    float alpha = 1;
    float beta = 0.f;

    gpuErrchk(cudaMemcpyAsync(d_sy, cudafphi_vars->h_y + batch_size*cudafphi_vars->stream_number*cudafphi_vars->n_subjects, sizeof(float)*cudafphi_vars->n_subjects*n_voxels, cudaMemcpyHostToDevice, stream));


    cudafphi_connectivity(cudafphi_vars->shared_device_vars.h2, cudafphi_vars->shared_device_vars.indicator, cudafphi_vars->shared_device_vars.pvals, d_sy, d_F,cudafphi_vars->shared_device_vars.hat,
    			cudafphi_vars->shared_device_vars.evectors, d_h2, d_indicator, d_pvals,  compute_F_vars,
    			cudafphi_vars->shared_device_vars.aux_vars, compute_h2_vars, pval_vars, cudafphi_vars->shared_device_vars.pmatrix,
    			cudafphi_vars->shared_device_vars.covariates,cudafphi_vars->shared_device_vars.get_pval, d_boolean_score,
    			n_voxels, cudafphi_vars->n_subjects, cudafphi_vars->n_permutations,
    			cudafphi_vars->stream_number, stream, handle);




    gpuErrchk(cudaStreamSynchronize(stream));
    gpuErrchk(cudaStreamDestroy(stream));

    cublasErrchk(cublasDestroy_v2(handle));



	gpuErrchk(cudaFree(compute_F_vars.d_Y));


	gpuErrchk(cudaFree(compute_F_vars.mean_or_sigma));



	gpuErrchk(cudaFree(compute_h2_vars.d_A));

	gpuErrchk(cudaFree(compute_h2_vars.d_B));

	gpuErrchk(cudaFree(compute_h2_vars.d_C));

	gpuErrchk(cudaFree(compute_h2_vars.d_D));

	gpuErrchk(cudaFree(compute_h2_vars.d_E));

	gpuErrchk(cudaFree(compute_h2_vars.d_Sigma_P));

	gpuErrchk(cudaFree(compute_h2_vars.d_score));

	gpuErrchk(cudaFree(compute_h2_vars.d_theta));

	gpuErrchk(cudaFree(compute_h2_vars.d_weights));


	gpuErrchk(cudaFree(compute_h2_vars.d_Sigma_E));

	gpuErrchk(cudaFree(compute_h2_vars.d_Sigma_A));

	if(cudafphi_vars->shared_device_vars.get_pval){
		gpuErrchk(cudaFree(pval_vars.syP));

		gpuErrchk(cudaFree(pval_vars.d_F));

		gpuErrchk(cudaFree( pval_vars.d_score));

		gpuErrchk(cudaFree( pval_vars.d_Ts));

		gpuErrchk(cudaFree( pval_vars.d_a));

		gpuErrchk(cudaFree( pval_vars.d_sigmaP));

	}


	gpuErrchk(cudaFree(d_sy));


	gpuErrchk(cudaFree(d_F));


	gpuErrchk(cudaFree(d_indicator));



	gpuErrchk(cudaFree(d_h2));


	gpuErrchk(cudaFree(d_boolean_score));


	if(cudafphi_vars->shared_device_vars.get_pval){
		gpuErrchk(cudaFree(d_pvals));
	}

    int * success = new int;
    *success = 1;
    return (void*)success;

}





