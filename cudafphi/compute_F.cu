#include "cudafPHI.cuh"
#include <stdio.h>
//#include <host_config.h>
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

static __global__ void calculate_means(float * d_sy, float * mean,  size_t n_voxels,size_t n_subjects){

  size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;
  size_t voxel = threadIdx.y + blockIdx.y*blockDim.y;

  __shared__ float shared_means[BLOCK_SIZE_1];
  size_t tIdx = threadIdx.x;
  if(rowIdx < n_subjects && voxel < n_voxels){
      shared_means[tIdx] = d_sy[voxel*n_subjects + rowIdx];
  }else{
      shared_means[tIdx] = 0.f;
  }

  __syncthreads();
  for(unsigned int stride = BLOCK_SIZE_1/2 ; stride > 0 ; stride >>=1){

	if(tIdx < stride && (rowIdx+stride < n_subjects)&& voxel < n_voxels){
	      shared_means[tIdx] += shared_means[tIdx + stride];
	}

	__syncthreads();
  }


  if((threadIdx.x == 0) && (rowIdx < n_subjects)&& voxel < n_voxels){
      atomicAdd(&mean[voxel], shared_means[0]);
  }
}


static __global__ void demean_columns(float * d_Y, float * d_sy, float * mean, size_t n_voxels, size_t n_subjects){

  size_t rowIdx = threadIdx.x +blockDim.x*blockIdx.x;
  size_t voxel = threadIdx.y + blockDim.y*blockIdx.y;


  if(rowIdx < n_subjects && voxel < n_voxels){
      d_Y[rowIdx + voxel*n_subjects] = d_sy[voxel*n_subjects + rowIdx] - mean[voxel]/float(n_subjects);
  }
}

static __global__ void calculate_sigma(float * d_sy, float * sigma,  size_t n_voxels,size_t n_subjects){

  size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;
  size_t voxel = threadIdx.y + blockIdx.y*blockDim.y;

  __shared__ float shared_means[BLOCK_SIZE_1];
  size_t tIdx = threadIdx.x;
  if(rowIdx < n_subjects && voxel < n_voxels){
      float value = d_sy[voxel*n_subjects + rowIdx];
      shared_means[tIdx] = value*value;
  }else{
      shared_means[tIdx] = 0.f;
  }

  __syncthreads();
  for(unsigned int stride = BLOCK_SIZE_1/2 ; stride > 0 ; stride >>=1){

	if(tIdx < stride && (rowIdx+stride < n_subjects)&& voxel < n_voxels){
	      shared_means[tIdx] += shared_means[tIdx + stride];
	}

	__syncthreads();
  }


  if((threadIdx.x == 0) && (rowIdx < n_subjects)&& voxel < n_voxels){
      atomicAdd(&sigma[voxel], shared_means[0]);
  }
}

static __global__ void standardize(float * d_Y,  float * sigma, size_t n_voxels, size_t n_subjects){

  size_t rowIdx = threadIdx.x +blockDim.x*blockIdx.x;
  size_t voxel = threadIdx.y + blockDim.y*blockIdx.y;


  if(rowIdx < n_subjects && voxel < n_voxels){
      d_Y[rowIdx + voxel*n_subjects] = d_Y[rowIdx + voxel*n_subjects]/sqrt(sigma[voxel]/float(n_subjects));
  }
}

static __global__ void calculate_SY(float * d_SY, const float * d_Y, const float * d_evectors, size_t n_voxels, size_t n_subjects){

  size_t rowIdx = blockDim.x*blockIdx.x + threadIdx.x;
  size_t colIdx = blockDim.y*blockIdx.y + threadIdx.y;

  size_t tIdx_x = threadIdx.x;
  size_t tIdx_y = threadIdx.y;

  float l_F = 0.f;


  __shared__  float shared_syP[BLOCK_SIZE_HAT][BLOCK_SIZE_HAT];
  __shared__  float shared_hat[BLOCK_SIZE_HAT][BLOCK_SIZE_HAT];
  for (size_t shared_kIdx = 0 ; shared_kIdx < gridDim.x; shared_kIdx++)
  {
      size_t kIdx_evec = shared_kIdx*BLOCK_SIZE_HAT + tIdx_y;
      size_t kIdx_SY = shared_kIdx*BLOCK_SIZE_HAT + tIdx_x;
      if(rowIdx >= n_subjects || kIdx_evec >= n_subjects){
          shared_hat[tIdx_x][tIdx_y] = 0.f;
      }else{
          shared_hat[tIdx_x][tIdx_y] = d_evectors[rowIdx*n_subjects + kIdx_evec];

      }

      if(colIdx >= n_voxels || kIdx_SY >= n_subjects){
          shared_syP[tIdx_x][tIdx_y] = 0.f;
      }else{
          shared_syP[tIdx_x][tIdx_y] = d_Y[kIdx_SY + colIdx*n_subjects];

      }

      __syncthreads();

      if((rowIdx < n_subjects) && (colIdx < (n_voxels))){
          for (int jIdx = 0; jIdx < BLOCK_SIZE_HAT; jIdx++){
   //   	local_hat_Idx = rowIdx*n_subjects + shared_kIdx*BLOCK_SIZE_HAT + k;
   //   	local_syP_Idx = colIdx*n_subjects + shared_kIdx*BLOCK_SIZE_HAT + k;
      	//if((local_hat_Idx < n_subjects*n_subjects) && (local_syP_Idx < n_subjects*(n_permutations + 1))){
      	l_F += shared_hat[tIdx_x][jIdx] * shared_syP[jIdx][tIdx_y];
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

/*
static __global__ void calculate_SY(float * d_SY, const float * d_Y, const float * d_evectors, size_t n_voxels, size_t n_subjects){
  size_t jIdx = threadIdx.x + blockIdx.x*blockDim.x;
  size_t rowIdx = threadIdx.y + blockIdx.y*blockDim.y;
  size_t voxel = blockIdx.z;
  __shared__ float shared_sy[BLOCK_SIZE_1];
  size_t tIdx = threadIdx.x;
  if(rowIdx < n_subjects  && jIdx < n_subjects){
      shared_sy[tIdx] = d_Y[jIdx*n_voxels + voxel]*d_evectors[rowIdx*n_subjects + jIdx];
  }else{
      shared_sy[tIdx] = 0.f;
  }

  __syncthreads();
  for(unsigned int stride = BLOCK_SIZE_1/2 ; stride > 0 ; stride >>=1){

	if(tIdx < stride && (jIdx+stride < n_subjects)){
	    shared_sy[tIdx] += shared_sy[tIdx + stride];
	}

	__syncthreads();
  }


  if((threadIdx.x == 0) && (jIdx < n_subjects)){
      atomicAdd(&d_SY[voxel*n_subjects + rowIdx], shared_sy[0]);
  }
}
*/

 /*
  * calcSXTSX

  * Input: d_sx (co-variate matrix), n_subjects(number of rows or subjects), n_covariates (number of co-variates or column_subjects in d_sx),
  * nD (number of data sets calculated)
  * Output: sxTsxI
  * After this function is called sxTsxI assumes the value of
  *  d_sx(tran_subjectsposed) * d_sx
  *  or written as
  *  d_sx' d_sx
  *  For further information about the linear algebra see computeCov function description
  */
static __global__ void calcSXTSX(float * d_sx, float * sxTsxI, size_t n_subjects,
		size_t n_covariates) {
	size_t rowIdx = threadIdx.x + blockIdx.x * blockDim.x;
	size_t colIdx1 = blockIdx.y;
	size_t colIdx2 = blockIdx.z;
	__shared__ float shared_sxTsxI[BLOCK_SIZE_1];


	if(rowIdx < n_subjects){
	    shared_sxTsxI[threadIdx.x] = d_sx[colIdx1*n_subjects + rowIdx]*d_sx[colIdx2*n_subjects + rowIdx];
	}else{
	    shared_sxTsxI[threadIdx.x] = 0.f;
	}
	__syncthreads();
	for(unsigned int stride = BLOCK_SIZE_1/2; stride > 0 ; stride >>= 1){
	    if((threadIdx.x < stride) && ((rowIdx + stride) < n_subjects)){
		shared_sxTsxI[threadIdx.x] += shared_sxTsxI[threadIdx.x + stride];
	    }
	    __syncthreads();
	}

	if((threadIdx.x == 0) && (rowIdx < n_subjects)){
	    atomicAdd(&sxTsxI[colIdx1*n_covariates + colIdx2], shared_sxTsxI[0]);
	}
}



/*
 * computeInverse
 * Input: A (n_subjects x n_subjects matrix stored on device), n_subjects (number of rows and column_subjects), and
 * nD (number of data sets in calculation)
 * Output: A (A^-1 overwrites the original values of A)
 *
 * The algorithm below relies on LU factorization to perform matrix inversion.  All calculation is
 * done on the device except kernels and memory tran_subjectsfers take place through loops on the CPU.
 *
 * As of 10/23/2014 original author has not tested whether this algorithm truly is faster than matrix inversion_subjects
 * using only the CPU.  But estimates that it should be faster for a large n_subjects and more so for a large nD.  Improvements
 * are always welcome and changes that would increase computation time should most likely start with using streams.  Also
 * speed and memory could be improved by making adjustments such as not allocating d_L as an n_subjects x n_subjects matrix but in_subjectstead a
 * matrix that contain_subjects only the lower elements since the upper elements are assumed to be 0.  With all that said this function is
 * called only once in the entire program (during calcRes) and the sizes of the matrices involved here would be n_covariates or the number of
 * column_subjects in the covariate matrix.
 *
 * d_U = A (Done through cudaMemcpy)
 * diagonal of d_L is set equal to 1 and all other entries are 0 GPU- InitialL
 *
 * d_L is gradually turned into a lower triangular matrix through a loop on CPU and a kernel
 * on the GPU.  GPU- setL
 *
 * d_U is gradually turned into a upper triangular matrix through the same loop on CPU calling kernel setL. GPU- setU
 *
 * d_M is gradually solved for from the equation: LM = I where I is the identity matrix.  M = U A^-1 GPU-solveL
 *
 * d_B is gradually solved for from the equation U d_B = M.  d_B = A^-1 and is copied over A taken as input. GPU- solveU
 *
 */
static int computeInverse(float * A, size_t n_subjects) {

  float * h_A = new float[n_subjects*n_subjects];
  char * functionName = "computeInverse";
  printError(cudaMemcpy(h_A, A, sizeof(float) * n_subjects * n_subjects,cudaMemcpyDeviceToHost),
             functionName,"cudaMemcpy-A to h_A", __LINE__);


/*
  if(n_subjects == 1){
      h_A[0] = 1.f/h_A[0];
      cudaMemcpy(A, h_A, sizeof(float)*n_subjects*n_subjects, cudaMemcpyHostToDevice);
      return 1;
  }else {
      float a = h_A[0];
      float b = h_A[1];
      float c = h_A[2];
      float d = h_A[3];

      float det = a*d - b*c;
      h_A[0] = d/det;
      h_A[1] = -c/det;
      h_A[2] = -b/det;
      h_A[3] = a/det;
      cudaMemcpy(A, h_A, sizeof(float)*n_subjects*n_subjects, cudaMemcpyHostToDevice);
      return 1;
  }*/
/*
  for(size_t j = 0; j < n_subjects ; j++){

      D[j] = h_A[j*n_subjects + j];
      float sum_ij = 0.f;

      for(size_t k = 0 ; k < j ; k++){
	  sum_ij += L[j*n_subjects + k]*L[j*n_subjects + k]*D[k];
      }
      D[j] = D[j] - sum_ij;
      for(size_t i = (j + 1) ; i < n_subjects ; i++){
	  float sum_ik = 0.f;
	  for(size_t k = 0 ; k < j ; k++){
	      sum_ik += L[k*n_subjects + i]*L[k*n_subjects + j]*D[k];
	  }
	  L[j*n_subjects + i] = (h_A[j*n_subjects + i] - sum_ik)/D[j];
      }
  }*/

  for(size_t i = 0 ; i < n_subjects ; i++){
      for(size_t j = 0 ; j < i; j++){
	  float sum = h_A[j*n_subjects + i];
	  for(size_t k = 0 ; k < j ;k++){
	      sum -= h_A[k*n_subjects + j]*h_A[j*n_subjects + k];
	  }
	  h_A[j*n_subjects + i] = sum/h_A[j*n_subjects + j];
      }
      for(size_t j = i ; j < n_subjects; j++){
	  float sum = h_A[j*n_subjects + i];
	  for(size_t k = 0 ; k < i ;k++){
	      sum -= h_A[k*n_subjects + j]*h_A[j*n_subjects + k];
	  }
	  h_A[j*n_subjects + i] = sum;
      }

  }


  float * M = new float[n_subjects*n_subjects];

  for(size_t j = 0 ; j < n_subjects ; j++){
      for(size_t i = 0 ; i < n_subjects ; i++){
	  double sum = 0.0;
	  if(i == j)
	    sum = 1.0;
	  for(size_t k = 0 ; k < i; k++){
	      sum -= M[j*n_subjects + k]*h_A[k*n_subjects + i];
	  }
	  M[j*n_subjects + i] = sum;
      }
  }




  float * h_AI = new float[n_subjects*n_subjects];
  for(size_t j = 0;  j < n_subjects; j++){
      for(size_t i = 0; i < n_subjects ; i++){
	  double sum = M[j*n_subjects + i];
	  for(size_t k = 0; k < i ; k++){
	      sum -= h_AI[j*n_subjects + k]*h_A[k*n_subjects + i];
	  }
	  h_AI[j*n_subjects + i] = sum/h_A[i*n_subjects + i];
      }
  }


/*
  float a = h_A[0];
  float b = h_A[1];
  float c = h_A[2];
  float d = h_A[3];

  float det = a*d - b*c;

  a = a/det;
  b = b/det;
  c = c/det;
  d = d/det;

  h_AI[0] = d;
  h_AI[1] = -c;
  h_AI[2] = -b;
  h_AI[3] = d;*
 // std::cout << a << " " << b << " " << c << " " << d << "\n";*/


  delete [] M;

  for(unsigned int i = 0; i < n_subjects ; i++){
      for(unsigned int j = 0 ; j < n_subjects ; j++){
	  printf("%f ", h_AI[j*n_subjects + i]);
      }
      printf("\n");
  }


  printError(cudaMemcpy(A, h_AI, sizeof(float)*n_subjects*n_subjects, cudaMemcpyHostToDevice),
             functionName, "cudaMemcpy-h_AI to A", __LINE__);

  delete [] h_AI;


  return 1;

}


/*
 * Input: d_sx (co-variate matrix), sxTsxI ( equals (d_sx'd_sx)^-1 ), n_subjects (number of subjects or rows),
 * n_covariates (number of co-variates or column_subjects of d_sx), nD (number of datasets calculated
 * Output: d_hat (hat matrix)
 * d_hat = I (Identity matrix) - d_sx * sxTsxI * d_sx(tran_subjectsposed)
 * For further information about the linear algebra see computeCov function description
 */
static __global__ void calcHat(float * d_sx, float * d_sxTsxI, float * d_hat,
		size_t n_subjects, size_t n_covariates) {

	size_t rowIdx = threadIdx.x + blockDim.x * blockIdx.x;
	size_t colIdx = blockIdx.y;



	if ((rowIdx < n_subjects) && (colIdx < n_subjects)) {
	    float value = 0.f;
	    for(size_t tIdx1 = 0 ; tIdx1 < n_covariates; tIdx1++){
		float local_trait = 0.f;
		for (size_t tIdx2 = 0; tIdx2 < n_covariates; tIdx2++) {
		    local_trait += d_sx[tIdx2 * n_subjects + rowIdx] * d_sxTsxI[tIdx1 * n_covariates + tIdx2];
		}
		value += local_trait*d_sx[tIdx1*n_subjects + colIdx];
	    }


	    if(rowIdx == colIdx){
		value = 1.f - value;
	    }else{
		value = -value;
	    }

	    d_hat[colIdx*n_subjects + rowIdx] = value;

	}

}


static __global__ void calculate_sum(float * d_B, float * d_Y, float * d_X,
                               size_t j, unsigned int rowIdx, unsigned int colIdx){


  size_t jIdx = threadIdx.x + blockIdx.x*blockDim.x;
  __shared__ float shared_sum[BLOCK_SIZE_1];

 // size_t j = gridDim.x;
 // size_t rowIdx = points.x;
  //size_t colIdx = points.y;

  if(jIdx < j){
     shared_sum[threadIdx.x] = d_Y[colIdx*j + jIdx] * d_X[rowIdx*j + jIdx];
 }else{
     shared_sum[threadIdx.x] = 0.f;
  }
  __syncthreads();
  for(unsigned int stride = BLOCK_SIZE_1/2; stride > 0 ; stride >>= 1){
      if((threadIdx.x < stride) && (jIdx + stride < j)){
	  shared_sum[threadIdx.x] += shared_sum[threadIdx.x + stride];
      }
      __syncthreads();
  }

  if((threadIdx.x == 0) &&(jIdx < j)){
      atomicAdd(&d_B[colIdx*j + rowIdx], shared_sum[0]);
  }


}
/*
static __global__ void calculate_F(float * d_SY, float * d_hat, float * d_F,
		size_t voxel, size_t n_subjects){
  size_t rowIdx = threadIdx.x + blockDim.x*blockIdx.x;
  __shared__ float shared_sy[BLOCK_SIZE_1];
  __shared__ float shared_hat[BLOCK_SIZE_1];
  size_t tIdx = threadIdx.x;
//  size_t tIdx_y = threadIdx.y;
  float l_F = 0.f;
  size_t l_syIdx;
  size_t l_hatIdx;
  for(size_t shared_jIdx = 0; shared_jIdx < gridDim.x; shared_jIdx++){
      l_syIdx = shared_jIdx*BLOCK_SIZE_1 + tIdx;
      l_hatIdx = rowIdx*n_subjects + shared_jIdx*BLOCK_SIZE_1 + tIdx;

      if(l_syIdx <  n_subjects){
	shared_sy[tIdx] = d_SY[voxel*n_subjects + l_syIdx];
      }else{
	shared_sy[tIdx] = 0.f;
      }

      if(l_hatIdx <  n_subjects*n_subjects){
	shared_sy[tIdx] = d_hat[l_hatIdx];
      }else{
	shared_sy[tIdx] = 0.f;
      }

      __syncthreads();

      if(rowIdx < n_subjects){
	  for(size_t jIdx = 0 ; jIdx < BLOCK_SIZE_1 ; jIdx++){
	      l_F += shared_sy[jIdx]*shared_hat[jIdx];
	  }
      }
      __syncthreads();
  }

  if(rowIdx >= n_subjects)
    return;

  d_F[voxel*n_subjects + rowIdx] = l_F*l_F;


}*/

static __global__ void calculate_F(const float * d_SY, const float * d_hat, float * d_F, size_t n_voxels,
                                   size_t n_subjects){

  size_t rowIdx = blockDim.x*blockIdx.x + threadIdx.x;
  size_t colIdx = blockDim.y*blockIdx.y + threadIdx.y;

  size_t tIdx_x = threadIdx.x;
  size_t tIdx_y = threadIdx.y;

  float l_F = 0.f;


  __shared__  float shared_syP[BLOCK_SIZE_HAT][BLOCK_SIZE_HAT];
  __shared__  float shared_hat[BLOCK_SIZE_HAT][BLOCK_SIZE_HAT];
  for (size_t shared_kIdx = 0 ; shared_kIdx < gridDim.x; shared_kIdx++)
    {
	 size_t kIdx_hat = shared_kIdx*BLOCK_SIZE_HAT + tIdx_y;
	 size_t kIdx_SY = shared_kIdx*BLOCK_SIZE_HAT + tIdx_x;
	// size_t local_hat_Idx = rowIdx*n_subjects + shared_kIdx*BLOCK_SIZE_HAT + tIdx_y;
	// size_t local_syP_Idx = colIdx*n_subjects + shared_kIdx*BLOCK_SIZE_HAT + tIdx_x;

      if(kIdx_hat >= n_subjects || rowIdx >= n_subjects){
          shared_hat[tIdx_x][tIdx_y] = 0.f;
      }else{
          shared_hat[tIdx_x][tIdx_y] = d_hat[rowIdx*n_subjects + kIdx_hat];

      }

      if(colIdx >= n_voxels || kIdx_SY >= n_subjects){
          shared_syP[tIdx_x][tIdx_y] = 0.f;
      }else{
          shared_syP[tIdx_x][tIdx_y] = d_SY[colIdx*n_subjects + kIdx_SY];

      }

      __syncthreads();

      if((rowIdx < n_subjects) && (colIdx < (n_voxels))){
          for (int jIdx = 0; jIdx < BLOCK_SIZE_HAT; jIdx++){
   //   	local_hat_Idx = rowIdx*n_subjects + shared_kIdx*BLOCK_SIZE_HAT + k;
   //   	local_syP_Idx = colIdx*n_subjects + shared_kIdx*BLOCK_SIZE_HAT + k;
      	//if((local_hat_Idx < n_subjects*n_subjects) && (local_syP_Idx < n_subjects*(n_permutations + 1))){
      	l_F += shared_hat[tIdx_x][jIdx] * shared_syP[jIdx][tIdx_y];
      	//}
          }
      }

      __syncthreads();

  }

  if((rowIdx >= n_subjects) || (colIdx >= (n_voxels))){
    return;
  }

  d_F[colIdx*n_subjects + rowIdx] = l_F;
}

/*
static __global__ void calculate_F(const float * d_SY, const float * d_hat, float * d_F,
		size_t n_voxels, size_t n_subjects) {

  size_t jIdx = threadIdx.x + blockIdx.x*blockDim.x;
  size_t rowIdx = threadIdx.y + blockIdx.y*blockDim.y;
  size_t voxel = blockIdx.z;
  __shared__ float shared_F[BLOCK_SIZE_1];
  size_t tIdx = threadIdx.x;
  if(rowIdx < n_subjects  && jIdx < n_subjects){
      shared_F[tIdx] = d_SY[jIdx + n_subjects*voxel]*d_hat[rowIdx*n_subjects + jIdx];
  }else{
      shared_F[tIdx] = 0.f;
  }

  __syncthreads();
  for(unsigned int stride = BLOCK_SIZE_1/2 ; stride > 0 ; stride >>=1){

	if(tIdx < stride && (jIdx+stride < n_subjects)){
	    shared_F[tIdx] += shared_F[tIdx + stride];
	}

	__syncthreads();
  }


  if((threadIdx.x == 0) && (jIdx < n_subjects)){
      atomicAdd(&d_F[voxel*n_subjects + rowIdx], shared_F[0]);
  }
}*/

static __global__ void set_F(float * d_F, float * d_sy, size_t n_voxels, size_t n_subjects){

  size_t rowIdx = blockIdx.x*blockDim.x + threadIdx.x;
  size_t voxel = blockIdx.y*blockDim.y + threadIdx.y;
  if(rowIdx < n_subjects){
      d_F[voxel*n_subjects + rowIdx] = d_sy[voxel*n_subjects + rowIdx];
  }


}

static __global__ void set_Zero(float * d_ptr, size_t n_subjects, size_t n_voxels){
  size_t rowIdx = blockIdx.x*blockDim.x + threadIdx.x;
  size_t voxel = blockIdx.y*blockDim.y + threadIdx.y;
  if(rowIdx < n_subjects){
      d_ptr[voxel*n_subjects + rowIdx] = 0.f;
  }
}

int compute_F(float * d_hat, float* d_sy, float *d_evectors, float * d_F,
              compute_F_variables vars, bool covariates, size_t n_subjects,
               size_t n_voxels, cudaStream_t stream) {


	const char * functionName = "compute_F.cu";

	int success;




	dim3 blockSize(BLOCK_SIZE_1, 1, 1);

	dim3 gridSize(ceil(float(n_subjects)/float(BLOCK_SIZE_1)), n_voxels, 1);

	dim3 gridSize2(ceil(float(n_subjects)/float(BLOCK_SIZE_1)), n_subjects, n_voxels);

	dim3 blockSizeHat(BLOCK_SIZE_HAT, BLOCK_SIZE_HAT, 1);
	dim3 gridSizeHat(ceil(float(n_subjects)/float(BLOCK_SIZE_HAT)), ceil(float(n_voxels)/float(BLOCK_SIZE_HAT)), 1);

	//for(unsigned int voxel = 0 ; voxel < n_voxels ; voxel++){

	    calculate_means<<<gridSize, blockSize, 0, stream>>>(d_sy, vars.mean, n_voxels, n_subjects);

	    demean_columns<<<gridSize, blockSize, 0, stream>>>(vars.d_Y, d_sy, vars.mean, n_voxels, n_subjects);
	    printError(cudaMemsetAsync(vars.mean, 0, sizeof(float)*n_voxels, stream ),
	               functionName, "cudaMemset-means", __LINE__);

	    calculate_sigma<<<gridSize, blockSize, 0, stream>>>(d_sy, vars.mean, n_voxels, n_subjects);

	    standardize<<<gridSize, blockSize, 0, stream>>>(vars.d_Y, vars.mean, n_voxels, n_subjects);
	  //  set_Zero<<<gridSize, blockSize, 0, stream>>>(d_sy, n_subjects, n_voxels);
	    calculate_SY<<<gridSizeHat, blockSizeHat, 0, stream>>>(d_sy, vars.d_Y, d_evectors, n_voxels, n_subjects);
	    if(covariates){
		calculate_F<<<gridSizeHat, blockSizeHat, 0, stream>>>(d_sy, d_hat, d_F, n_voxels, n_subjects);
	    }else{
		set_F<<<gridSize, blockSize, 0, stream>>>(d_F, d_sy, n_voxels, n_subjects);
	    }

	//}





	return 1;
}
