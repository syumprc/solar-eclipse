#include "cudafphi.cuh"
#include <stdio.h>
#include <math_functions.h>
#include <iostream>

static __global__ void calculate_means(const float * d_sy, float * mean, const size_t n_voxels, const size_t n_subjects){


	/*cublasHandle_t handle;
	cublasCreate_v2(&handle);

	size_t vIdx = threadIdx.x;

	cublasSetPointerMode_v2(handle, CUBLAS_POINTER_MODE_DEVICE);

    cublasDasum_v2(handle, n_subjects, &d_sy[vIdx*n_subjects], 1, &mean[vIdx]);*/


		unsigned int tIdx = threadIdx.x;

		unsigned int rowIdx = blockIdx.x*blockDim.x + threadIdx.x;
		unsigned int vIdx = blockIdx.y*blockDim.y + threadIdx.y;


		extern __shared__ float shared_mean[];

		if(rowIdx < n_subjects && vIdx < n_voxels)
			shared_mean[tIdx] = d_sy[vIdx*n_subjects + rowIdx];
		else
			shared_mean[tIdx] = 0.0;




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
	  const float l_mean = mean[voxel]/float(n_subjects);
	  const float value  =  d_sy[voxel*n_subjects + rowIdx] - l_mean;
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

static __global__ void calculate_SY(float * d_SY, const float * d_Y, const float * d_evectors, const size_t n_voxels, const size_t n_subjects){

  const size_t rowIdx = blockDim.x*blockIdx.x + threadIdx.x;
  const size_t colIdx = blockDim.y*blockIdx.y + threadIdx.y;

  const size_t tIdx_x = threadIdx.x;
  const size_t tIdx_y = threadIdx.y;

  float l_F = 0.f;


  extern __shared__  float shared_data[];
  //__shared__  float shared_hat[BLOCK_SIZE_HAT][BLOCK_SIZE_HAT];
  bool should_return_1 = false;
  bool should_return_2 = false;
  for (size_t shared_kIdx = 0 ; shared_kIdx < gridDim.x; shared_kIdx++)
  {
      size_t kIdx_evec = shared_kIdx*blockDim.x + tIdx_y;
      size_t kIdx_SY = shared_kIdx*blockDim.x + tIdx_x;
      if(rowIdx >= n_subjects || kIdx_evec >= n_subjects){
          //shared_hat[tIdx_x][tIdx_y] = 0.f;
    	  shared_data[tIdx_x*blockDim.x + tIdx_y] = 0.f;
    	  should_return_1 = true;
      }else{
          //shared_hat[tIdx_x][tIdx_y] = d_evectors[rowIdx*n_subjects + kIdx_evec];
    	  shared_data[tIdx_x*blockDim.x + tIdx_y] =  d_evectors[rowIdx*n_subjects + kIdx_evec];
      }

      if(colIdx >= n_voxels || kIdx_SY >= n_subjects){
    	  shared_data[blockDim.x*blockDim.x + tIdx_x*blockDim.x + tIdx_y] = 0.f;
    	  should_return_2 = true;
      }else{
    	  shared_data[blockDim.x*blockDim.x + tIdx_x*blockDim.x + tIdx_y] = d_Y[kIdx_SY + colIdx*n_subjects];

      }

    ///  if(should_return_1 && should_return_2)
   // 	  return;


      __syncthreads();

      if((rowIdx < n_subjects) && (colIdx < (n_voxels))){
          for (int jIdx = 0; jIdx < blockDim.x; jIdx++){
   //   	local_hat_Idx = rowIdx*n_subjects + shared_kIdx*BLOCK_SIZE_HAT + k;
   //   	local_syP_Idx = colIdx*n_subjects + shared_kIdx*BLOCK_SIZE_HAT + k;
      	//if((local_hat_Idx < n_subjects*n_subjects) && (local_syP_Idx < n_subjects*(n_permutations + 1))){
      	l_F +=  shared_data[tIdx_x*blockDim.x + jIdx] * shared_data[blockDim.x*blockDim.x + jIdx*blockDim.x + tIdx_y];
      	//}
          }
      }

      __syncthreads();

  }

  if((rowIdx >= n_subjects) || (colIdx >= (n_voxels))){
    return;
  }

  d_SY[colIdx*n_subjects + rowIdx] = l_F;
}



int compute_F(const float * d_hat, float* d_sy, const float *d_evectors, float * d_F,
              compute_F_variables vars, bool covariates, size_t n_subjects,
               size_t n_voxels, cudaStream_t stream, cublasHandle_t handle) {





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

	dim3 blockSize_mult(16,16,1);
	dim3 gridSize_mult(ceil(float(n_subjects)/float(16.0)), ceil(float(n_voxels)/float(16.0)), 1);
    gpuErrchk(cudaMemsetAsync(vars.mean_or_sigma, 0, sizeof(float)*n_voxels, stream ));





    calculate_means<<<gridSize_mean, blockSize_mean, blockSize_n_subjects*sizeof(float), stream>>>(d_sy, vars.mean_or_sigma, n_voxels, n_subjects);
    gpuErrchk(cudaPeekAtLastError());
	demean_columns<<<gridSize_mean, blockSize_mean, 0, stream>>>(vars.d_Y, d_sy, vars.mean_or_sigma, n_voxels, n_subjects);
    gpuErrchk(cudaPeekAtLastError());

	float alpha = 1.0;
	float beta = 0.0;
   gpuErrchk(cudaMemsetAsync(vars.mean_or_sigma, 0, sizeof(float)*n_voxels, stream ));

	calculate_sigma<<<gridSize_mean, blockSize_mean, sizeof(float)*1024, stream>>>(vars.d_Y, vars.mean_or_sigma, n_voxels, n_subjects);
	gpuErrchk(cudaPeekAtLastError());
	calculate_inverse_normal<<<gridSize_mean, blockSize_mean, 0, stream>>>(vars.d_Y, vars.mean_or_sigma, n_voxels, n_subjects);
	gpuErrchk(cudaPeekAtLastError());

//	gpuErrchk(cudaStreamSynchronize(stream));
   cublasErrchk(cublasSgemm_v2(handle, CUBLAS_OP_T, CUBLAS_OP_N, n_subjects,n_voxels, n_subjects, &alpha, d_evectors, n_subjects, vars.d_Y, n_subjects, &beta, d_sy, n_subjects));

//	gpuErrchk(cudaStreamSynchronize(stream));

//	calculate_SY<<<gridSize_mult, blockSize_mult, sizeof(float)*16*16*2, stream>>>(d_sy,  vars.d_Y,  d_evectors, n_voxels, n_subjects);

//	gpuErrchk(cudaStreamSynchronize(stream));

   if(covariates){

	    cublasErrchk(cublasSgemm_v2(handle, CUBLAS_OP_N, CUBLAS_OP_N, n_subjects,n_voxels, n_subjects, &alpha, d_hat, n_subjects, d_sy, n_subjects, &beta, d_F, n_subjects));

   }else{
	   gpuErrchk(cudaMemcpyAsync(d_F, d_sy, sizeof(float)*n_subjects*n_voxels, cudaMemcpyDeviceToDevice, stream));

   }


   return 1;
}
