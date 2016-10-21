#include "cudafPHI.cuh"
#include "cuda_runtime.h"
#include <stdio.h>





//surface<void, 2> surfD;
static inline void printError(cudaError_t error, const char * functionMess, const char * opMess, int Line) {
	if (error != 0){
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

__global__ void hat_syP( float* d_F, const float* d_hat, const float* d_syP,
                         size_t n_subjects, size_t n_permutations, size_t voxel , bool * d_boolean_score)
{

    if(d_boolean_score[voxel] == false){
	return;
    }
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
	 size_t kIdx_syP = shared_kIdx*BLOCK_SIZE_HAT + tIdx_x;

        if((rowIdx >= n_subjects || kIdx_hat >= n_subjects)){
            shared_hat[tIdx_x][tIdx_y] = 0.f;
        }else{
            shared_hat[tIdx_x][tIdx_y] = d_hat[rowIdx*n_subjects + kIdx_hat];

        }

        if(kIdx_syP >= n_subjects || colIdx >= (n_permutations + 1)){
            shared_syP[tIdx_x][tIdx_y] = 0.f;
        }else{
            shared_syP[tIdx_x][tIdx_y] = d_syP[colIdx*n_subjects + kIdx_syP];

        }

        __syncthreads();

        if((rowIdx < n_subjects) && (colIdx < (n_permutations + 1))){
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

    if((rowIdx >= n_subjects) || (colIdx >= (n_permutations + 1))){
      return;
    }

    d_F[colIdx*n_subjects + rowIdx] = l_F*l_F;

}




/*
 * Input: d_res (residual vector), d_cov (regressed sy of covariate matrix from ordinary least squares see computeCov for
 * further information), d_pmatrix (matrix of random indices), is_covariate (determines if d_a covariate matrix was provided),
 * n_subjects (number of subjects or rows), nD (number of data sets handled), and n_permutations (number of permutations)
 * Output: syP (permutation matrix)
 */

static __global__ void set_syP(float * syP, const float * d_res, const float * d_sy,
                               const unsigned int * d_indices, unsigned int voxel,
 const size_t n_permutations, const size_t n_subjects, size_t n_voxels, bool * d_boolean_score){
  if(d_boolean_score[voxel] == false)
    return;
  const unsigned int row_Idx = threadIdx.x + blockDim.x*blockIdx.x;
  const unsigned int p_Idx = blockIdx.y;

  if((p_Idx < (n_permutations + 1)) && (row_Idx < n_subjects)){
      if(p_Idx == 0){
	  float value = d_sy[row_Idx +  n_subjects*voxel];
	  syP[row_Idx] = value;
	  //surf2Dwrite(value, surfD, row_Idx*sizeof(float), 0);
      }else{
	  size_t Idx = d_indices[(p_Idx - 1)*n_subjects + row_Idx];
	  float value = d_res[voxel*n_subjects + Idx] - (d_res[voxel*n_subjects + row_Idx] - d_sy[voxel*n_subjects + row_Idx]);

	 // surf2Dwrite(value, surfD, row_Idx*sizeof(float), p_Idx);
	  syP[p_Idx*n_subjects + row_Idx] = value;
      }

  }

}

static __global__ void calculate_sum(float * B, const float * A, const float * X, size_t rows, uint3 previous_threadIdx,
                                     uint3 previous_blockIdx, dim3 previous_blockDim){

  size_t jIdx = threadIdx.x + blockDim.x*blockIdx.x;
  if(jIdx >= rows)
    return;

  size_t colIdx = previous_blockIdx.x*previous_blockDim.x + previous_threadIdx.x;
  size_t rowIdx = previous_blockIdx.y*previous_blockDim.y + previous_threadIdx.y;
  size_t tIdx = threadIdx.x;
  __shared__ float shared_sum[BLOCK_SIZE_1];

  shared_sum[tIdx] = A[colIdx*rows + jIdx]*X[rowIdx*rows + jIdx];
  __syncthreads();

  for(unsigned int stride = BLOCK_SIZE_1/2 ; stride > 0 ; stride >>= 1){
      if((tIdx < stride) && ((jIdx + stride) < rows))
	shared_sum[tIdx] += shared_sum[tIdx + stride];

      __syncthreads();
  }

  if(tIdx == 0)
    atomicAdd(B, shared_sum[0]);


}




/*static __global__ void hat_syP(float * d_F, const float * d_syP, const float * d_hat,
                               size_t n_subjects, size_t n_permutations, size_t voxel, bool * d_boolean_score) {
  if(d_boolean_score[voxel] == false)
    return;

	const size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;
	const size_t colIdx = blockIdx.y;


	__shared__ float shared_P[16][16];
	int i = threadIdx.x;
	int j = threadIdx.y;

	shared_P[i][j] = 0.f;
	if((rowIdx >= n_subjects) || (colIdx >= (n_permutations + 1)))
	    return;

	const size_t syP_place = colIdx*n_subjects;
	const size_t d_hat_place = rowIdx*n_subjects;
	//const float * l_hat = &d_hat[rowIdx*n_subjects];
	for(int w = 0 ; w < ceil(float(n_subjects)/16.f) ; w++){
	    size_t jIdx = 16*w + j;
	    if(jIdx < n_subjects){
		shared_P[i][j] += (d_syP[syP_place + jIdx])*(d_hat[d_hat_place + jIdx]);
	    }

	}
	__syncthreads();

	  for(int stride = 16/2 ; stride > 0 ; stride >>= 1){
	      if((j < stride))
		shared_P[i][j] += shared_P[i][j + stride];

	      __syncthreads();
	  }
	  if(j > 0)
	    return;
	  float value = shared_P[i][0];
	  d_F[colIdx*n_subjects + rowIdx] = value*value;
	//  atomicAdd(&d_F[colIdx*n_subjects + rowIdx], shared_P[i][0]);



	//__shared__ float shared_F[4];
	//int tIdx = threadIdx.y;
	//shared_F[tIdx] = d_syP[colIdx*n_subjects + jIdx]*d_hat[rowIdx*n_subjects + jIdx];
//	cudaStream_t s;
//	cudaStreamCreateWithFlags(&s, cudaStreamNonBlocking);
//	calculate_sum<<<ceil(float(n_subjects)/float(BLOCK_SIZE_1)), BLOCK_SIZE_1, 0, s>>>(&d_F[colIdx*n_subjects + rowIdx], d_syP, d_hat, n_subjects, threadIdx, blockIdx, blockDim);



	//size_t rowIdx = blockIdx.z*blockDim.z + threadIdx.z;
	//size_t jtIdx = threadIdx.x;

//	float sum = 0.f;
//	for(size_t jIdx = 0; jIdx < n_subjects ; jIdx++){
	    //float value;
	    //surf2Dread(&value, surfD, jIdx*sizeof(float), colIdx);
//	    sum += d_syP[colIdx*n_subjects + jIdx]*d_hat[rowIdx*n_subjects + jIdx];
//	}
//	d_F[colIdx*n_subjects + rowIdx] = sum*sum;
	//size_t coltIdx = threadIdx.y;
	/*__shared__ float shared_syP[32][16];
	__shared__ float shared_hat[16][16];
	float value = 0.f;

	float * local_F = new float[32*16];



	for(size_t r = 0 ; r < 16 ; r++){
	    for(size_t c = 0 ; c < 32 ; c ++){
		local_F[c*16 + r] = d_F[]
	    }
	}

	for(size_t m = 0 ; m < ceil(n_subjects/16.f) ; m++){

	    if ((jIdx < n_subjects) && (rowIdx < n_subjects) && (colIdx < (n_permutations + 1))) {
		   // float value;
		   // surf2Dread(&value,  surfD, jIdx*sizeof(float), colIdx, cudaBoundaryModeTrap);
		shared_syP[coltIdx][jtIdx] = d_syP[colIdx*n_subjects + jIdx];
		shared_hat[coltIdx][jtIdx] = d_hat[rowIdx * n_subjects + jIdx];
	    }else{
		shared_syP[coltIdx][jtIdx] = 0.f;
		shared_hat[coltIdx][jtIdx] = 0.f;
	    }

	    __syncthreads();




	}*/
/*	__syncthreads();
	for(unsigned int stride = 4/2 ; stride > 0 ; stride >>= 1){
	    if((threadIdx.x < stride) && ((jIdx + stride) < n_subjects))
	      shared_F[tIdx] += shared_F[tIdx + stride];

	    __syncthreads();
	}




	if(tIdx == 0){
	    atomicAdd(&d_F[colIdx*n_subjects + rowIdx], shared_F[0]);
	}

	//size_t pIdx = threadIdx.x + blockDim.x*blockIdx.x;
	//size_t rowIdx = blockIdx.y;
	    const float * l_syP = &d_syP[colIdx*n_subjects];
	    const float * l_hat = &d_hat[rowIdx*n_subjects];
	    float sum = 0.f;
	    for(unsigned int jIdx = 0 ; jIdx < n_subjects ; jIdx++){
		sum += (*l_syP++)*(*l_hat++);
	    }
	    d_F[colIdx*n_subjects + rowIdx] = sum*sum;




}*/

static __global__ void square_F(float * d_F, size_t n_subjects, size_t voxel, bool * d_boolean_score){
  if(d_boolean_score[voxel] == false)
    return;
  size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;
  size_t colIdx = blockIdx.y;

  if(rowIdx < (n_subjects)){
      float value = d_F[colIdx*n_subjects + rowIdx];
      d_F[colIdx*n_subjects + rowIdx] = value*value;
  }
}

static __global__ void set_F(float * d_F, const float * syP, size_t n_subjects, size_t n_permutations,size_t voxel, bool * d_boolean_score){
  if(d_boolean_score[voxel] == false)
    return;
  size_t pIdx = threadIdx.x + blockIdx.x*blockDim.x;
  size_t rowIdx = blockIdx.y;

  if(pIdx < (n_permutations + 1)){
      float value = syP[pIdx*n_subjects + rowIdx];
      d_F[pIdx*n_subjects + rowIdx] = value*value;
  }
}



static __global__ void calculate_sigmaP(const float * syP, float * d_sigmaP, size_t n_subjects,
		 size_t n_permutations, size_t voxel, bool * d_boolean_score) {
         if(d_boolean_score[voxel] == false)
           return;
	const size_t rowIdx = threadIdx.x + blockDim.x * blockIdx.x;
	size_t pIdx = blockIdx.y*blockDim.y + threadIdx.y;
	__shared__ float shared_sum[BLOCK_SIZE_1];
	if ((rowIdx < n_subjects) && (pIdx < n_permutations + 1)) {

	    shared_sum[threadIdx.x] = syP[pIdx*n_subjects + rowIdx];

	}else{

	    shared_sum[threadIdx.x] = 0.f;
	    return;
	}

	__syncthreads();

	for(unsigned int stride = BLOCK_SIZE_1/2 ; stride > 0 ; stride >>= 1){
	    if((threadIdx.x < stride) && ((rowIdx + stride) < n_subjects))
	      shared_sum[threadIdx.x] += shared_sum[threadIdx.x + stride];

	    __syncthreads();
	}


	if((threadIdx.x == 0) && (rowIdx < n_subjects)){
	    atomicAdd(&d_sigmaP[pIdx], shared_sum[0]);
	}
}

static __global__ void update_syP(float * syP,const float * d_sigmaP, size_t n_subjects,
		size_t n_permutations, size_t voxel,bool * d_boolean_score) {
        if(d_boolean_score[voxel] == false)
          return;
	size_t rowIdx = threadIdx.x + blockDim.x * blockIdx.x;
	size_t pIdx = blockIdx.y*blockDim.y + threadIdx.y;
	if ((rowIdx < n_subjects) && (pIdx < n_permutations + 1)) {
		float div_num = d_sigmaP[pIdx] / float(n_subjects);
		syP[pIdx*n_subjects +rowIdx] = ((syP[pIdx*n_subjects +rowIdx] / (div_num)));
	}
}



static __global__ void calculate_score(const float * d_Z, const float * syP, const float * d_sigmaP,
		float * d_score, size_t n_subjects, size_t n_permutations,size_t voxel, bool * d_boolean_score) {
        if(d_boolean_score[voxel] == false)
          return;
	size_t rowIdx = threadIdx.x + blockDim.x * blockIdx.x;
	size_t  pIdx = blockIdx.y*blockDim.y + threadIdx.y;
	size_t tIdx = threadIdx.x;
	__shared__ float shared_score[BLOCK_SIZE_1];

	if ((rowIdx < n_subjects) && (pIdx < (n_permutations+1))) {
		float div_num = d_sigmaP[pIdx]  / float(n_subjects);
		shared_score[tIdx] = (0.5 * d_Z[n_subjects + rowIdx] * (syP[pIdx*n_subjects + rowIdx]) - 1.f);

	}else{
	    shared_score[tIdx] = 0.f;
	    return;
	}
	__syncthreads();
	for(unsigned int stride = BLOCK_SIZE_1/2 ; stride > 0 ; stride >>= 1){

	    if((tIdx < stride) && ((rowIdx + stride) < n_subjects)){
	      shared_score[tIdx] += shared_score[tIdx + stride];
	    }

	    __syncthreads();
	}


	if((threadIdx.x == 0) && (rowIdx < n_subjects)){
	    atomicAdd(&d_score[pIdx], shared_score[0]);
	}
}


static __global__ void calculate_A(const float * syP, const float * d_Z,  float * d_a, size_t n_subjects,
		size_t n_permutations,size_t voxel, bool * d_boolean_score) {
         if(d_boolean_score[voxel] == false)
           return;
	size_t rowIdx = threadIdx.x + blockDim.x * blockIdx.x;
	size_t pIdx = blockIdx.y*blockDim.y + threadIdx.y;

	__shared__ float shared_a_0[BLOCK_SIZE_1];
	__shared__ float shared_a_1[BLOCK_SIZE_1];

	size_t tIdx = threadIdx.x;
	if ((rowIdx < n_subjects) && (pIdx < n_permutations +1)) {
		shared_a_0[tIdx] = syP[pIdx*n_subjects + rowIdx];
		shared_a_1[tIdx] = syP[pIdx*n_subjects + rowIdx] * d_Z[rowIdx + n_subjects];
	}else{
	    shared_a_0[tIdx] = 0.f;
	    shared_a_1[tIdx] = 0.f;
	    return;
	}

	__syncthreads();

	for(unsigned int stride = BLOCK_SIZE_1/2 ; stride > 0 ; stride >>= 1){

	    if((threadIdx.x < stride) && ((rowIdx + stride) < n_subjects )){
	      shared_a_0[tIdx] += shared_a_0[tIdx + stride];
	      shared_a_1[tIdx] += shared_a_1[tIdx + stride];
	    }

	    __syncthreads();
	}


	if((tIdx == 0) && (rowIdx < n_subjects)){
	    atomicAdd(&d_a[pIdx*2 ], shared_a_0[0]);
	    atomicAdd(&d_a[pIdx*2 + 1], shared_a_1[0]);
	}


}

/*
 * Input: d_a ((n_permutations + 1) x 2 matrix used to calculate d_Ts), d_b (d_b = (Z' Z)^-1 ) ,nD (number of data sets handled),
 * and n_permutations (number of permutations)
 * Output: d_Ts (Test statistic vector)
 *
 * c = d_a * d_b
 * d_Ts = 0.5 (d_a[column 0]*c[column 0] + d_a[column 1]*c[column 1])
 */
static __global__ void calculate_Ts(float * d_a, const float * d_sigmaP, const float * d_b, float * d_Ts,
                       size_t n_permutations, size_t n_subjects, size_t voxel, bool * d_boolean_score) {
	if(d_boolean_score[voxel] == false)
	  return;
	 size_t pIdx = threadIdx.x + blockDim.x * blockIdx.x;
	if ((pIdx < (n_permutations + 1))) {
	    float l_a_1 = d_a[pIdx*2];
	    float l_a_2 = d_a[pIdx*2 + 1];
/*
	   if(l_a_1 < 0.f)
	     l_a_1 = 0.f;

	    if(l_a_2 < 0.f)
	      l_a_2 = 0.f;*/


	    float row_0 = l_a_1*d_b[0] + l_a_2*d_b[1];


	    row_0 = row_0*l_a_1;
	    //row_0 = row_0*row_0;

	    float row_1 = l_a_1*d_b[2] + l_a_2*d_b[3];

	    row_1 = row_1*l_a_2;

	    d_Ts[pIdx] = 0.5*(row_1 + row_0);
	}
}


static __global__ void compare_score(const float * d_score, float * d_Ts, size_t n_permutations, size_t voxel, bool * d_boolean_score) {

        if(d_boolean_score[voxel] == false)
          return;
	size_t rowIdx = threadIdx.x + blockDim.x * blockIdx.x;
	if (rowIdx < (n_permutations + 1)) {
	    d_Ts[rowIdx] = d_score[rowIdx];
//	    if (d_score[rowIdx] < 0.f) {
//		d_Ts[rowIdx] = d_score[rowIdx];
//	    }
	}


}




/*
 * Input: d_Ts (Test statistic vector), nD (number of data sets handled), and
 * n_permutations (number of permutations)
 * Output: d_pval (pval statistic)
 *
 * if(d_Ts[i] >= T[0])
 *  d_pval += 1
 */
static __global__ void sum_Ts(const float * d_Ts, float * d_pval, size_t vIdx, size_t n_permutations, size_t voxel, bool * d_boolean_score){

  if(d_boolean_score[voxel] == false)
    return;

  size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;
//  size_t vIdx = blockIdx.y;
  size_t tIdx = threadIdx.x;

  __shared__ float shared_pval[BLOCK_SIZE_2];

  if (rowIdx < (n_permutations + 1)) {
      float comp_val = d_Ts[0];
      if(d_Ts[rowIdx] >= comp_val){
	  shared_pval[tIdx] = 1.f;
      }else{
	  shared_pval[tIdx] = 0.f;
      }
  }else{
      shared_pval[tIdx] = 0.f;
      return;
  }

  __syncthreads();

  for(size_t stride = BLOCK_SIZE_2/2 ; stride > 0 ; stride >>= 1){

      if(tIdx < stride  && ((rowIdx + stride) < (n_permutations + 1))){
	  shared_pval[tIdx] += shared_pval[tIdx + stride];
      }

      __syncthreads();
  }

  if((threadIdx.x == 0) && (rowIdx < (n_permutations + 1))){
      atomicAdd(&d_pval[vIdx], shared_pval[0]);
  }


}

/*
 * Input: d_pval (currently sum of one's that depend on the d_Ts vector), nD (number of data sets), and n_permutations
 * (number of permutations)
 * Output: d_pval (average the sum calculated in sumTs)
 */
static __global__ void mean_Ts(float * d_pval, float n_permutations, size_t n_voxels, bool * d_boolean_score) {



  size_t vIdx = blockIdx.x*blockDim.x + threadIdx.x;
  if(vIdx >= n_voxels)
    return;
  if(d_boolean_score[vIdx] == false){
      d_pval[vIdx] = 0.5f;
      return;
  }
  d_pval[vIdx] = d_pval[vIdx]/(n_permutations + 1.f);


}



/*
 * calcPval
 * Input: res (residual vector), cov (cov = sy - res), d_Z (auxiliary matrix), hat (hat matrix
 * hat = (I - sx (sx' sx)^-1 sx')), h_pmatrix (matrix of random indices), is_covariate (determines whether
 * or not d_a covariate matrix was provided), n_subjects (number of subjects or rows), nD (number of data sets handled),
 * n_permutations (number of permutations), and tileNum (tile number out of total data set)
 * Output: pval (pval associated with h^2 and sigma statistics)
 *
 * syP[column 0] = res       Sets first column equal to res----------------------------------------
 * syP[column i] = res[h_pmatrix[column i - 1]]  ((if co_variate) + cov) Set all other columns----GPU - setsyP
 * equal to d_a randomized version of res and if is_covariate = true then cov is added.--------------
 *
 * if(is_covariate = true)
 * syP = hat*syP - GPU- HatsyP
 *
 * syP = syP*.syP (*. is element wise multiplication) GPU- sqsyP
 *
 * d_b = Z' Z GPU- calcZTZ
 *
 * d_b = d_b^-1 GPU- ZTZInv
 *
 * d_sigmaP[i] = sum(syP[column i]) GPU- calcsigmaP
 *
 * syP = (syP/(d_sigmaP/n_subjects) - 1) GPU- updatesyP
 *
 * d_score[i] = 0.5 * d_Z[eigencolumn] * syP[column i]/(d_sigmaP[i]/n_subjects) GPU- calcScore
 *
 *
 *
 * d_a = d_Z' * syP GPU- calcA
 *
 * c =  d_a*d_b-------------------------------------------------------GPU- calcTs
 * d_Ts = 0.5*(c[column 0]*.d_a[column 0]  + c[column 1] *. d_a[column 1]) --/
 *
 * if(d_score[i] <= 0)-----/
 *  d_Ts[i] = 0 -----------GPU- compareScore
 *
 *  pval = 1-----------------/
 *  if(d_Ts[i] >= d_Ts[0])---GPU- sumTs
 *   pval += 1 --------------/
 *
 *  pval = pval/(n_permutations + 1) GPU- meanTs
 */


int compute_pvals(float * d_sy, float * d_res, float * d_hat, float * d_sigma_E,
                  float * d_sigma_A, float * d_pvals, unsigned int * d_pmatrix,
                  bool * d_boolean_score, aux_variables aux_vars, pval_variables pval_vars, bool covariate,
                  size_t n_subjects, size_t n_voxels, size_t n_permutations, cudaStream_t stream) {


	const char * functionName = "compute_pvals.cu";
	size_t freeMem;
	size_t totalMem;

	//cudaChannelFormatDesc channelDesc =   cudaCreateChannelDesc(8, 8, 8, 8,
	  //                                                          cudaChannelFormatKindFloat);

	   //Alternatively
	//cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

	//cudaArray *d_arr;   cudaMallocArray(&d_arr, &channelDesc, n_subjects, (n_permutations + 1), cudaArraySurfaceLoadStore);
	//printf("here\n");

	///printError(cudaBindSurfaceToArray(surfD, d_arr), functionName, "cudaBindSurfaceToArray", __LINE__);

	//printf("here\n");





	size_t Idx;
	dim3 blockSize1(BLOCK_SIZE_2, 1, 1);
	dim3 blockSize2(BLOCK_SIZE_1, 1, 1);
	dim3 gridSize1(ceil(float(n_permutations + 1)/float(BLOCK_SIZE_2)), n_subjects, 1);

	dim3 gridSize2(ceil(float(n_subjects)/float(BLOCK_SIZE_1)), n_permutations + 1, 1);

	int gridSizeHatx;
	int blockSizeHatx;

	 dim3 blockSizeHat(BLOCK_SIZE_HAT , BLOCK_SIZE_HAT, 1);

	 dim3 gridSizeHat((n_subjects + (BLOCK_SIZE_HAT - 1))/BLOCK_SIZE_HAT ,(n_permutations + 1 + (BLOCK_SIZE_HAT - 1))/BLOCK_SIZE_HAT, 1);


	for(size_t voxel = 0 ; voxel < n_voxels ; voxel++){

	 //   if(h_boolean_score[voxel]){


		set_syP<<<gridSize2, blockSize2, 0, stream>>>(pval_vars.syP, (const float *)d_res, (const float *)d_sy, (const unsigned int *)d_pmatrix, voxel, n_permutations, n_subjects , n_voxels, d_boolean_score);

		//printError(cudaMemsetAsync(pval_vars.d_F, 0, sizeof(float) * (n_permutations + 1) * n_subjects, stream),
		//functionName, "cudaMemset-d_sigmaP", __LINE__);
		if(covariate){
		//    hat_syP<<<gridSizeHat, blockSizeHat, 0, stream>>>(pval_vars.d_F, (const float *)pval_vars.syP, (const float *)d_hat, n_subjects, n_permutations, voxel, d_boolean_score);
		    hat_syP<<<gridSizeHat, blockSizeHat, 0, stream>>>( pval_vars.d_F,  d_hat, pval_vars.syP,  n_subjects,  n_permutations, voxel, d_boolean_score);
		   // square_F <<<gridSize2, blockSize2, 0 , stream>>>(pval_vars.d_F, n_subjects);
		}else{
		    set_F<<<gridSize1, blockSize1, 0, stream>>>(pval_vars.d_F,(const float*)pval_vars.syP, n_subjects, n_permutations, voxel, d_boolean_score);
		}



		printError(cudaMemsetAsync(pval_vars.d_sigmaP, 0, sizeof(float) * (n_permutations + 1), stream),
			functionName, "cudaMemset-d_sigmaP", __LINE__);

		calculate_sigmaP<<<gridSize2, blockSize2, 0, stream>>>((const float*)pval_vars.d_F, pval_vars.d_sigmaP, n_subjects, n_permutations, voxel, d_boolean_score);

		update_syP<<<gridSize2, blockSize2, 0, stream>>>(pval_vars.d_F,(const float*) pval_vars.d_sigmaP, n_subjects, n_permutations, voxel,d_boolean_score);

		printError(cudaMemsetAsync(pval_vars.d_score, 0, sizeof(float) * (n_permutations + 1), stream),
			functionName, "cudaMemset-d_score", __LINE__);

		calculate_score<<<gridSize2, blockSize2, 0, stream>>>((const float*)aux_vars.d_Z,(const float*) pval_vars.d_F, (const float *)pval_vars.d_sigmaP,
		    pval_vars.d_score, n_subjects, n_permutations, voxel, d_boolean_score);


	//	printError(cudaMemsetAsync(pval_vars.d_a, 0, sizeof(float) * (n_permutations + 1) * 2, stream) ,
	//		functionName, "cudaMemset-d_a", __LINE__);

	//	calculate_A<<<gridSize2, blockSize2, 0, stream>>>((const float*)pval_vars.d_F, (const float*)aux_vars.d_Z, pval_vars.d_a, n_subjects, n_permutations,voxel, d_boolean_score);

	//	printError(cudaMemsetAsync(pval_vars.d_Ts, 0.f, sizeof(float) * (n_permutations + 1), stream),
	//		functionName, "cudaMemset-d_Ts", __LINE__);

	//	calculate_Ts<<<ceil(float(n_permutations + 1)/float(BLOCK_SIZE_2)), BLOCK_SIZE_2, 0, stream>>>(pval_vars.d_a, pval_vars.d_sigmaP,
	//	    (const float*) aux_vars.d_ZTZI, pval_vars.d_Ts, n_permutations, n_subjects, voxel, d_boolean_score);

		compare_score<<<ceil(float(n_permutations + 1)/float(BLOCK_SIZE_2)), BLOCK_SIZE_2, 0, stream>>>((const float*)pval_vars.d_score, pval_vars.d_Ts, n_permutations, voxel,d_boolean_score);

		sum_Ts<<<ceil(float(n_permutations + 1)/float(BLOCK_SIZE_2)), BLOCK_SIZE_2, 0, stream>>>((const float *)pval_vars.d_Ts, d_pvals, voxel, n_permutations, voxel,d_boolean_score);


	 //   }


	}
	mean_Ts<<<ceil(float(n_voxels)/float(BLOCK_SIZE_3)), BLOCK_SIZE_3, 0, stream >>>(d_pvals, n_permutations, n_voxels, d_boolean_score);









	return 1;

}
