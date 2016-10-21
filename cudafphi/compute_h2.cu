#include "cudafPHI.cuh"
#include <stdio.h>
#include <cuda_runtime.h>


static inline void printError(cudaError_t error, const char * functionMess, const char * opMess, int Line) {
	if (error == 0){
		return;
	}else{
		//fclose(stderr);
		//freopen("cudafPHI.err", "w", stderr);
		printf("%s %s Line - %i : %s \n",functionMess , opMess, Line ,cudaGetErrorString(error) );
		//std::cerr << functionMess << " " << opMess << " Line- " << Line << ": "
		//		<< cudaGetErrorString(error) << "\n";
		//fclose(stderr);
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


static __global__ void calculate_d_theta(float * d_theta, float * d_F, float * d_Z, float * d_ZTZI, size_t voxel, size_t n_subjects){

  size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;
  size_t tIdx = threadIdx.x;
  __shared__ float shared_theta_0[BLOCK_SIZE_1];
  __shared__ float shared_theta_1[BLOCK_SIZE_1];

  if (rowIdx < n_subjects) {

      float d_Z_0 = d_Z[rowIdx];
      float d_Z_1 = d_Z[rowIdx + n_subjects];
      float F = d_F[rowIdx + voxel*n_subjects];
      F = F*F;
      shared_theta_0[tIdx] = F*(d_Z_0*d_ZTZI[0] + d_ZTZI[2]*d_Z_1);
      shared_theta_1[tIdx] = F*(d_Z_0*d_ZTZI[1] + d_ZTZI[3]*d_Z_1);

    //  printf("%f \n",shared_theta_0[threadIdx.x]);
  }else{
      shared_theta_0[tIdx] = 0.f;
      shared_theta_1[tIdx] = 0.f;
  }
  __syncthreads();
  for(unsigned int stride = BLOCK_SIZE_1/2 ; stride > 0 ; stride >>=1){

	if(tIdx < stride && (rowIdx+stride < n_subjects)){
	      shared_theta_0[tIdx] += shared_theta_0[tIdx + stride];
	      shared_theta_1[tIdx] += shared_theta_1[tIdx + stride];
	}

	__syncthreads();
  }


  if((tIdx == 0) && (rowIdx < n_subjects)){
      atomicAdd(&d_theta[0], shared_theta_0[0]);
      atomicAdd(&d_theta[1], shared_theta_1[0]);
  }

}


static __global__ void calculate_weights(float * d_weights, float * d_Z, float * d_theta, size_t n_subjects){


  size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;

  if(rowIdx < n_subjects){
      float theta_0 = d_theta[0];
      float theta_1 = d_theta[1];
     if(theta_0 < 0.f)
  	theta_0 = 0.f;

     if(theta_1 < 0.f)
  	theta_1 = 0.f;

      float weight = theta_0 + d_Z[rowIdx + n_subjects]*theta_1;
      d_weights[rowIdx] = 1.f/(weight*weight);

  }


}



static __global__ void calculate_sum_A_D_Sigma_P(float * d_A, float * d_D, float * d_Sigma_P, float * d_weights,
                                          float * d_F, size_t voxel, size_t n_subjects){

  size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;
  size_t tIdx = threadIdx.x;
  __shared__ float shared_A[BLOCK_SIZE_1];
  __shared__ float shared_D[BLOCK_SIZE_1];
  __shared__ float shared_Sigma_P[BLOCK_SIZE_1];

  if(rowIdx < n_subjects){
      float weight = d_weights[rowIdx];
      float F = d_F[voxel*n_subjects + rowIdx];
      F = F*F;
      shared_A[tIdx] = weight;
      shared_D[tIdx] = weight*F;
      shared_Sigma_P[tIdx] = F;
  }else{
      shared_A[tIdx] = 0.f;
      shared_D[tIdx] = 0.f;
      shared_Sigma_P[tIdx] = 0.f;
  }

  __syncthreads();
  for(unsigned int stride = BLOCK_SIZE_1/2 ; stride > 0 ; stride >>=1){

	if(tIdx < stride && (rowIdx+stride < n_subjects)){
	      shared_A[tIdx] += shared_A[tIdx + stride];
	      shared_D[tIdx] += shared_D[tIdx + stride];
	      shared_Sigma_P[tIdx] += shared_Sigma_P[tIdx + stride];
	}

	__syncthreads();
  }


  if((tIdx == 0) && (rowIdx < n_subjects)){
      atomicAdd(d_A, shared_A[0]);
      atomicAdd(d_D, shared_D[0]);
      atomicAdd(d_Sigma_P, shared_Sigma_P[0]);
  }
}

static __global__ void calculate_sum_B_C_E(float * d_B, float * d_C, float * d_E, float * d_weights, float * d_F,
                                    float * d_Z, size_t voxel, size_t n_subjects){

  size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;
  size_t tIdx = threadIdx.x;
  __shared__ float shared_B[BLOCK_SIZE_1];
  __shared__ float shared_C[BLOCK_SIZE_1];
  __shared__ float shared_E[BLOCK_SIZE_1];

  if(rowIdx < n_subjects){
      float lambda = d_Z[rowIdx + n_subjects];
      float F = d_F[voxel*n_subjects + rowIdx];
      F = F*F;
      float weight = d_weights[rowIdx];
      shared_B[tIdx] = weight*lambda;
      shared_C[tIdx] = weight*lambda*lambda;
      shared_E[tIdx] = F*weight*lambda;
  }else{
      shared_B[tIdx] = 0.f;
      shared_C[tIdx] = 0.f;
      shared_E[tIdx] = 0.f;
  }

  __syncthreads();
  for(unsigned int stride = BLOCK_SIZE_1/2 ; stride > 0 ; stride >>=1){

	if(tIdx < stride && (rowIdx+stride < n_subjects)){
	      shared_B[tIdx] += shared_B[tIdx + stride];
	      shared_C[tIdx] += shared_C[tIdx + stride];
	      shared_E[tIdx] += shared_E[tIdx + stride];
	}

	__syncthreads();
  }


  if((threadIdx.x == 0) && (rowIdx < n_subjects)){
      atomicAdd(d_B, shared_B[0]);
      atomicAdd(d_C, shared_C[0]);
      atomicAdd(d_E, shared_E[0]);
  }
}

static __global__ void calculate_score(float * d_score, float * d_Z, float * d_F, float * d_Sigma_P, size_t voxel, size_t n_subjects){

  size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;

  size_t tIdx = threadIdx.x;

  __shared__ float shared_score[BLOCK_SIZE_1];

  if(rowIdx < n_subjects){
      float F = d_F[voxel*n_subjects + rowIdx];
      F = F*F;
      shared_score[tIdx] = 0.5f*d_Z[rowIdx + n_subjects]*(F/(*d_Sigma_P/float(n_subjects)) - 1.f);
  }else{
      shared_score[tIdx] = 0.f;
  }

  __syncthreads();

  for(unsigned int stride = BLOCK_SIZE_1/2 ; stride > 0 ; stride >>=1){

	if(tIdx < stride && (rowIdx+stride < n_subjects)){
	      shared_score[tIdx] += shared_score[tIdx + stride];
	}

	__syncthreads();
  }
  if((tIdx == 0) && (rowIdx < n_subjects)){
      atomicAdd(d_score, shared_score[0]);
  }
}



static __global__ void calculate_h2(float * d_h2, float * d_indicator, float * d_Sigma_P, float * d_score, bool * d_boolean_score, float * d_A, float * d_B, float * d_C, float * d_D, float * d_E, size_t voxel){

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
	  if(h2r <= 0.f){
	      h2r = 0.f;
	      d_boolean_score[voxel] = false;
	  }else{
	      d_boolean_score[voxel] = true;
	  }
	  d_h2[voxel] = h2r;

	  *d_Sigma_P = sigma_E + sigma_A;

      }
}

/*
static __global__ void calculate_weights_SE(float * d_weights, bool * d_boolean_score, float * d_Z, float * d_h2, float * d_Sigma_P, size_t n_subjects, size_t voxel){

  if(d_boolean_score[voxel] == false){
      return;
  }
  size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;

  if(rowIdx < n_subjects){
      float Sigma_A = d_h2[voxel]*(*d_Sigma_P);
      float Sigma_E = *d_Sigma_P - Sigma_A;
      if(Sigma_A < 0.f)
	Sigma_A = 0.f;

      if(Sigma_E < 0.f)
	Sigma_E = 0.f;

      float weight = Sigma_E + d_Z[rowIdx + n_subjects]*Sigma_A;
      d_weights[rowIdx] = 1.f/(weight*weight);

  }


}


static __global__ void calculate_sum_A_B_C(float * d_A, float * d_B, float * d_C, float * d_weights,
                                           float * d_Z, bool * d_boolean_score, size_t n_subjects, size_t voxel){

  if(d_boolean_score[voxel] == false)
    return;


  size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;
  size_t tIdx = threadIdx.x;
  __shared__ float shared_A[BLOCK_SIZE_1];
  __shared__ float shared_B[BLOCK_SIZE_1];
  __shared__ float shared_C[BLOCK_SIZE_1];

  if(rowIdx < n_subjects){
      float weight = d_weights[rowIdx];
      float lambda = d_Z[rowIdx + n_subjects];
      shared_A[tIdx] = weight;
      shared_B[tIdx] = weight*F;
      shared_C[tIdx] = F;
  }else{
      shared_A[tIdx] = 0.f;
      shared_B[tIdx] = 0.f;
      shared_C[tIdx] = 0.f;
  }

  __syncthreads();
  for(unsigned int stride = BLOCK_SIZE_1/2 ; stride > 0 ; stride >>=1){

	if(tIdx < stride && (rowIdx+stride < n_subjects)){
	    shared_A[tIdx] += shared_A[tIdx + stride];
	    shared_B[tIdx] += shared_B[tIdx + stride];
	    shared_C[tIdx] += shared_C[tIdx + stride];
	}

	__syncthreads();
  }


  if((tIdx == 0) && (rowIdx < n_subjects)){
      atomicAdd(d_A, shared_A[0]);
      atomicAdd(d_B, shared_B[0]);
      atomicAdd(d_C, shared_C[0]);
  }


}
*/
int compute_h2(float * d_F, float * d_h2,  float * d_indicator, bool * d_boolean_score,
               compute_h2_variables vars, aux_variables aux_vars,
               size_t n_subjects, size_t n_voxels, cudaStream_t stream) {


	const char * functionName = "compute_h2.cu";



	dim3 blockSize(BLOCK_SIZE_1, 1, 1);
	dim3 gridSize(ceil(float(n_subjects) / float(BLOCK_SIZE_1)), 1, 1);



	blockSize.x = BLOCK_SIZE_1;
	gridSize.x = ceil(float(n_subjects)/float(BLOCK_SIZE_1));
	blockSize.y = 1;
	gridSize.y = 1;

	for(size_t voxel = 0; voxel < n_voxels ; voxel++){


	    printError(cudaMemsetAsync(vars.d_theta, 0, sizeof(float)*2, stream),
			functionName, "cudaMemset-d_theta", __LINE__);


	    calculate_d_theta<<<gridSize, blockSize, 0, stream>>>(vars.d_theta, d_F, aux_vars.d_Z, aux_vars.d_ZTZI, voxel, n_subjects);


	    calculate_weights<<<gridSize, blockSize, 0, stream>>>(vars.d_weights, aux_vars.d_Z, vars.d_theta, n_subjects);

	    printError(cudaMemsetAsync(vars.d_A, 0, sizeof(float), stream),
			functionName, "cudaMemset-d_A", __LINE__);

	    printError(cudaMemsetAsync(vars.d_D, 0, sizeof(float), stream),
			functionName, "cudaMemset-d_D", __LINE__);


	    printError(cudaMemsetAsync(vars.d_Sigma_P, 0, sizeof(float), stream),
			functionName, "cudaMemset-d_Sigma_P", __LINE__);


	    calculate_sum_A_D_Sigma_P<<<gridSize, blockSize, 0, stream>>>(vars.d_A, vars.d_D, vars.d_Sigma_P,
		vars.d_weights, d_F, voxel, n_subjects);

	    printError(cudaMemsetAsync(vars.d_B, 0, sizeof(float), stream),
			functionName, "cudaMemset-d_B", __LINE__);

	    printError(cudaMemsetAsync(vars.d_C, 0, sizeof(float), stream),
			functionName, "cudaMemset-d_C", __LINE__);

	    printError(cudaMemsetAsync(vars.d_E, 0, sizeof(float),stream),
			functionName, "cudaMemset-d_E", __LINE__);

	    calculate_sum_B_C_E<<<gridSize, blockSize, 0 ,stream>>>(vars.d_B, vars.d_C, vars.d_E, vars.d_weights, d_F, aux_vars.d_Z, voxel, n_subjects);


	    printError(cudaMemsetAsync(vars.d_score, 0, sizeof(float), stream),
			functionName, "cudaMemset-d_score", __LINE__);

	    calculate_score<<<gridSize, blockSize, 0, stream>>>(vars.d_score, aux_vars.d_Z, d_F, vars.d_Sigma_P, voxel, n_subjects);



	    calculate_h2<<<1, 1, 0, stream>>>(d_h2, d_indicator, vars.d_Sigma_P, vars.d_score, d_boolean_score, vars.d_A, vars.d_B, vars.d_C, vars.d_D, vars.d_E, voxel);




	}




	return 1;
}
