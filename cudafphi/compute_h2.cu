#include "cudafphi.cuh"
#include <stdio.h>
#include <cuda_runtime.h>


static __global__ void calculate_d_theta(float * d_theta, float * d_F, const float * d_Z, float * d_ZTZI, size_t voxel, size_t n_subjects){

  size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;
  size_t tIdx = threadIdx.x;
  extern __shared__ float shared_theta[];


  if (rowIdx < n_subjects) {

      float d_Z_0 = d_Z[rowIdx];
      float d_Z_1 = d_Z[rowIdx + n_subjects];
      float F = d_F[rowIdx + voxel*n_subjects];
      F = F*F;

      shared_theta[tIdx*2] = F*(d_Z_0*d_ZTZI[0] + d_ZTZI[2]*d_Z_1);
      shared_theta[tIdx*2 + 1] = F*(d_Z_0*d_ZTZI[1] + d_ZTZI[3]*d_Z_1);

    //  printf("%f \n",shared_theta_0[threadIdx.x]);
  }else{
	  shared_theta[tIdx*2] = 0.f;
	  shared_theta[tIdx*2 + 1] = 0.f;
  }
  __syncthreads();
  for(unsigned int stride = blockDim.x/2 ; stride > 0 ; stride >>=1){

	if(tIdx < stride && (rowIdx+stride < n_subjects)){
		shared_theta[tIdx*2] += shared_theta[(tIdx + stride)*2];
		shared_theta[tIdx*2 + 1] += shared_theta[(tIdx + stride)*2 + 1];
	}

	__syncthreads();
  }


  if((tIdx == 0) && (rowIdx < n_subjects)){
      atomicAdd(&d_theta[0], shared_theta[0]);
      atomicAdd(&d_theta[1], shared_theta[1]);
  }



}


static __global__ void calculate_weights(float * d_weights, const float * d_Z, float * d_theta, size_t n_subjects){


  size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;

  if(rowIdx < n_subjects){
      float theta_0 = d_theta[0];
      float theta_1 = d_theta[1];
      if(theta_0 < 0.f){
    	  theta_0 = 0.f;

      }

      if(theta_1 < 0.f){
    	  theta_1 = 0.f;

      }

      float weight = d_Z[rowIdx]*theta_0 + d_Z[rowIdx + n_subjects]*theta_1;

      if(weight == 0.f){
    	  d_weights[rowIdx] = 0.f;
      }else{
    	  d_weights[rowIdx] = 1.f/(weight*weight);
      }




  }


}



static __global__ void calculate_sum_A_D_Sigma_P(float * d_A, float * d_D, float * d_Sigma_P, float * d_weights,
                                          float * d_F, size_t voxel, size_t n_subjects){

  size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;
  size_t tIdx = threadIdx.x;

  extern __shared__ float shared_data[];

  if(rowIdx < n_subjects){
      float weight = d_weights[rowIdx];
      float F = d_F[voxel*n_subjects + rowIdx];
      F = F*F;

      shared_data[tIdx*3] = weight;
      shared_data[tIdx*3 + 1] = weight*F;
      shared_data[tIdx*3 + 2] = F;

  }else{
      shared_data[tIdx*3] = 0.f;
      shared_data[tIdx*3 + 1] = 0.f;
      shared_data[tIdx*3 + 2] = 0.f;

  }

  __syncthreads();
  for(unsigned int stride = blockDim.x/2; stride  > 0 ; stride >>=1){

	  if(tIdx < stride && (rowIdx +stride < n_subjects)){
	      shared_data[tIdx*3] += shared_data[(tIdx + stride)*3];
	      shared_data[tIdx*3 + 1] += shared_data[(tIdx + stride)*3 + 1];
	      shared_data[tIdx*3 + 2] += shared_data[(tIdx + stride)*3 + 2];
	  }

	  __syncthreads();
  }


  if((tIdx == 0) && (rowIdx < n_subjects)){
      atomicAdd(d_A, shared_data[0]);
      atomicAdd(d_D, shared_data[1]);
      atomicAdd(d_Sigma_P, shared_data[2]);
  }
}

static __global__ void calculate_sum_B_C_E(float * d_B, float * d_C, float * d_E, float * d_weights, float * d_F,
                                    const float * d_Z, size_t voxel, size_t n_subjects){

  size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;
  size_t tIdx = threadIdx.x;
  extern __shared__ float shared_data[];

  if(rowIdx < n_subjects){

      float lambda = d_Z[rowIdx + n_subjects];
      float F = d_F[voxel*n_subjects + rowIdx];
      F = F*F;
      float weight = d_weights[rowIdx];

      shared_data[tIdx*3] = weight*lambda;
      shared_data[tIdx*3 + 1] = weight*lambda*lambda;
      shared_data[tIdx*3 + 2] = F*weight*lambda;

  }else{
      shared_data[tIdx*3] = 0.f;
      shared_data[tIdx*3 + 1] = 0.f;
      shared_data[tIdx*3 + 2] = 0.f;

  }

  __syncthreads();
  for(unsigned int stride = blockDim.x/2; stride  > 0 ; stride >>=1){

	  if(tIdx < stride && (rowIdx+stride < n_subjects)){
	      shared_data[tIdx*3] += shared_data[(tIdx + stride)*3];
	      shared_data[tIdx*3 + 1] += shared_data[(tIdx + stride)*3 + 1];
	      shared_data[tIdx*3 + 2] += shared_data[(tIdx + stride)*3 + 2];
	  }

	  __syncthreads();
  }


  if((threadIdx.x == 0) && (rowIdx < n_subjects)){
      atomicAdd(d_B, shared_data[0]);
      atomicAdd(d_C, shared_data[1]);
      atomicAdd(d_E, shared_data[2]);
  }
}

static __global__ void calculate_score(float * d_score, const float * d_Z, float * d_F, float * d_Sigma_P, size_t voxel, size_t n_subjects){

  size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;

  size_t tIdx = threadIdx.x;

  extern __shared__ float shared_score[];

  if(rowIdx < n_subjects){
      float F = d_F[voxel*n_subjects + rowIdx];
      F = F*F;
      shared_score[tIdx] = 0.5f*d_Z[rowIdx + n_subjects]*(F/(*d_Sigma_P/float(n_subjects)) - 1.f);
  }else{
      shared_score[tIdx] = 0.f;
  }

  __syncthreads();

  for(unsigned int stride = blockDim.x/2 ; stride > 0 ; stride >>=1){

	if(tIdx < stride && (rowIdx+stride < n_subjects)){
	      shared_score[tIdx] += shared_score[tIdx + stride];
	}

	__syncthreads();
  }
  if((tIdx == 0) && (rowIdx < n_subjects)){
      atomicAdd(d_score, shared_score[0]);
  }
}



static __global__ void calculate_h2(float * d_h2, float * d_indicator, float * d_Sigma_A, float * d_Sigma_E, float * d_score, bool * d_boolean_score, float * d_A, float * d_B, float * d_C, float * d_D, float * d_E, size_t voxel){

      d_indicator[voxel] = *d_score;
      float A = *d_A;
      float B = *d_B;
      float C = *d_C;
      float D = *d_D;
      float E = *d_E;

      float sigma_E = (C*D - B*E)/(A*C - B*B);
      float sigma_A = (A*E - B*D)/(A*C - B*B);

      if(sigma_E < 0.f)
    	  sigma_E = 0.f;

      if(sigma_A < 0.f)
    	  sigma_A = 0.f;

      float h2r;

      if((sigma_E + sigma_A) == 0.f){
    	  d_h2[voxel] = 0.f;
	  	  d_boolean_score[voxel] = false;
      }else{

    	  if(*d_score < 0){
    		  h2r = 0.f;
    	  }
    	  h2r = sigma_A/(sigma_E + sigma_A);
    	  if(h2r <= 0.f || h2r != h2r){
    		  h2r = 0.f;
	      	  d_boolean_score[voxel] = false;
    	  }else{
    		  d_boolean_score[voxel] = true;
    	  }
	  	  d_h2[voxel] = h2r;

	  	  d_Sigma_A[voxel] = sigma_A;
	  	  d_Sigma_E[voxel] = sigma_E;

      }
}





int compute_h2(float * d_F, float * d_h2,  float * d_indicator, bool * d_boolean_score,
               compute_h2_variables vars, aux_variables aux_vars,
               size_t n_subjects, size_t n_voxels, cudaStream_t stream) {


	const char * functionName = "compute_h2.cu";

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
	dim3 gridSize(ceil(float(n_subjects)/float(blockSize_n_subjects)), 1, 1);



	for(size_t voxel = 0; voxel < n_voxels ; voxel++){
	    gpuErrchk(cudaMemsetAsync(vars.d_theta, 0, sizeof(float)*2, stream));


	    calculate_d_theta<<<gridSize, blockSize, sizeof(float)*blockSize_n_subjects*2, stream>>>(vars.d_theta, d_F, aux_vars.d_Z, aux_vars.d_ZTZI, voxel, n_subjects);
	    gpuErrchk(cudaPeekAtLastError());

	    calculate_weights<<<gridSize, blockSize, 0, stream>>>(vars.d_weights, aux_vars.d_Z, vars.d_theta, n_subjects);
	    gpuErrchk(cudaPeekAtLastError());

	    gpuErrchk(cudaMemsetAsync(vars.d_A, 0, sizeof(float), stream));

	    gpuErrchk(cudaMemsetAsync(vars.d_D, 0, sizeof(float), stream));


	    gpuErrchk(cudaMemsetAsync(vars.d_Sigma_P, 0, sizeof(float), stream));


	    calculate_sum_A_D_Sigma_P<<<gridSize, blockSize, sizeof(float)*blockSize_n_subjects*3, stream>>>(vars.d_A, vars.d_D, vars.d_Sigma_P,
		vars.d_weights, d_F, voxel, n_subjects);
	    gpuErrchk(cudaPeekAtLastError());

	    gpuErrchk(cudaMemsetAsync(vars.d_B, 0, sizeof(float), stream));

	    gpuErrchk(cudaMemsetAsync(vars.d_C, 0, sizeof(float), stream));

	    gpuErrchk(cudaMemsetAsync(vars.d_E, 0, sizeof(float),stream));

	    calculate_sum_B_C_E<<<gridSize, blockSize, sizeof(float)*blockSize_n_subjects*3 ,stream>>>(vars.d_B, vars.d_C, vars.d_E, vars.d_weights, d_F, aux_vars.d_Z, voxel, n_subjects);
	    gpuErrchk(cudaPeekAtLastError());


	    gpuErrchk(cudaMemsetAsync(vars.d_score, 0, sizeof(float), stream));

	    calculate_score<<<gridSize, blockSize, sizeof(float)*blockSize_n_subjects, stream>>>(vars.d_score, aux_vars.d_Z, d_F, vars.d_Sigma_P, voxel, n_subjects);
	    gpuErrchk(cudaPeekAtLastError());



	    calculate_h2<<<1, 1, 0, stream>>>(d_h2, d_indicator, vars.d_Sigma_A, vars.d_Sigma_E, vars.d_score, d_boolean_score, vars.d_A, vars.d_B, vars.d_C, vars.d_D, vars.d_E, voxel);
	    gpuErrchk(cudaPeekAtLastError());




	}




	return 1;
}
