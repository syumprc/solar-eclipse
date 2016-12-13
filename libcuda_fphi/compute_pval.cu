#include "cudafphi.cuh"
#include "cuda_runtime.h"
#include <stdio.h>



__global__ void hat_syP( float* d_F, const float* d_hat, const float* d_syP,
                         size_t n_subjects, size_t n_permutations, size_t voxel , bool * d_boolean_score)
{

    if(d_boolean_score[voxel] == false){
	return;
    }

    const size_t rowIdx = blockDim.x*blockIdx.x + threadIdx.x;
    const size_t colIdx = blockDim.y*blockIdx.y + threadIdx.y;

    const size_t tIdx_x = threadIdx.x;
    const size_t tIdx_y = threadIdx.y;

    float l_F = 0.f;



#pragma unroll
    for (size_t shared_kIdx = 0 ; shared_kIdx < gridDim.x; shared_kIdx++)
      {

        __shared__  float shared_syP[BLOCK_SIZE_HAT][BLOCK_SIZE_HAT];
        __shared__  float shared_hat[BLOCK_SIZE_HAT][BLOCK_SIZE_HAT];
    	const size_t kIdx_hat = shared_kIdx*BLOCK_SIZE_HAT + tIdx_y;
    	const size_t kIdx_syP = shared_kIdx*BLOCK_SIZE_HAT + tIdx_x;
        if((rowIdx >= n_subjects)  || kIdx_hat >= n_subjects){
            shared_hat[tIdx_x][tIdx_y] = 0.f;
        }else{
            shared_hat[tIdx_x][tIdx_y] = d_hat[kIdx_hat*n_subjects + rowIdx];

        }

        if(colIdx >= (n_permutations + 1) || kIdx_syP >= n_subjects){
            shared_syP[tIdx_x][tIdx_y] = 0.f;
        }else{
            shared_syP[tIdx_x][tIdx_y] = d_syP[colIdx*n_subjects + kIdx_syP];

        }


        if(colIdx < (n_permutations + 1) && (rowIdx < n_subjects)){
			#pragma unroll (16)
            for (int jIdx = 0; jIdx < blockDim.x; jIdx++){
            	l_F += shared_hat[tIdx_x][jIdx] * shared_syP[jIdx][tIdx_y];

            }

        }

        __syncthreads();

    }


    d_F[colIdx*n_subjects + rowIdx] = l_F*l_F;

}




static __global__ void set_syP(float * syP, const float * d_res, const float * d_sy, const unsigned int * d_indices,
		unsigned int voxel, const size_t n_permutations, const size_t n_subjects, size_t n_voxels, bool * d_boolean_score){
  if(d_boolean_score[voxel] == false)
    return;

  const unsigned int row_Idx = threadIdx.x + blockDim.x*blockIdx.x;
  const unsigned int p_Idx = blockIdx.y;

  if((p_Idx < (n_permutations + 1)) && (row_Idx < n_subjects)){
      if(p_Idx == 0){
    	  const float value = d_sy[row_Idx +  n_subjects*voxel];
    	  syP[row_Idx] = value;
      }else{
    	  size_t Idx = d_indices[(p_Idx - 1)*n_subjects + row_Idx];

    	  const float value = d_res[Idx +  n_subjects*voxel] + (d_sy[voxel*n_subjects + row_Idx] - d_res[voxel*n_subjects + row_Idx]);
    	  syP[p_Idx*n_subjects + row_Idx] = value;
      }

  }

}






static __global__ void set_F(float * d_F, const float * syP, const unsigned int * pmatrix, size_t n_subjects, size_t n_permutations,size_t voxel, bool * d_boolean_score){
  if(d_boolean_score[voxel] == false)
    return;
  const size_t pIdx = threadIdx.y + blockIdx.y*blockDim.y;
  const size_t rowIdx = blockIdx.x*blockDim.x + threadIdx.x;

  if(pIdx < (n_permutations + 1)  && rowIdx < n_subjects){
	  const float value = syP[pIdx*n_subjects + rowIdx];
      d_F[pIdx*n_subjects + rowIdx] = value*value;
  }
}



static __global__ void calculate_sigmaP(const float * syP, float * d_sigmaP, size_t n_subjects,
		 size_t n_permutations, size_t voxel, bool * d_boolean_score) {
         if(d_boolean_score[voxel] == false)
           return;
	const size_t rowIdx = threadIdx.x + blockDim.x * blockIdx.x;
	const size_t pIdx = blockIdx.y*blockDim.y + threadIdx.y;
	extern __shared__ float shared_sum[];
	if ((rowIdx < n_subjects) && (pIdx < n_permutations + 1)) {

	    shared_sum[threadIdx.x] = syP[pIdx*n_subjects + rowIdx];

	}else{

	    shared_sum[threadIdx.x] = 0.f;
	    return;
	}

	__syncthreads();

	for(unsigned int stride = blockDim.x/2 ; stride > 0 ; stride >>= 1){
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
	const size_t rowIdx = threadIdx.x + blockDim.x * blockIdx.x;
	const size_t pIdx = blockIdx.y*blockDim.y + threadIdx.y;
	if ((rowIdx < n_subjects) && (pIdx < n_permutations + 1)) {
		const float div_num = d_sigmaP[pIdx] / float(n_subjects);
		const float val = (syP[pIdx*n_subjects +rowIdx]/div_num - 1.f);
		syP[pIdx*n_subjects +rowIdx]  = val;
	}
}



static __global__ void calculate_score(const float * d_Z, const float * syP, const float * d_sigmaP,
		float * d_score, size_t n_subjects, size_t n_permutations,size_t voxel, bool * d_boolean_score) {
        if(d_boolean_score[voxel] == false)
          return;
	size_t rowIdx = threadIdx.x + blockDim.x * blockIdx.x;
	size_t  pIdx = blockIdx.y*blockDim.y + threadIdx.y;
	size_t tIdx = threadIdx.x;
	extern __shared__ float shared_score[];

	if ((rowIdx < n_subjects) && (pIdx < (n_permutations+1))) {
		const float div_num = d_sigmaP[pIdx]  / float(n_subjects);
		const float score_val = d_Z[n_subjects + rowIdx];
		shared_score[tIdx] = (0.5f * score_val  * ( (syP[pIdx*n_subjects + rowIdx])/div_num - 1.f));

	}else{
	    shared_score[tIdx] = 0.f;
	    return;
	}
	__syncthreads();
	for(unsigned int stride = blockDim.x/2 ; stride > 0 ; stride >>= 1){

	    if((tIdx < stride) && ((rowIdx + stride) < n_subjects)){
	      shared_score[tIdx] += shared_score[tIdx + stride];
	    }

	    __syncthreads();
	}


	if((threadIdx.x == 0) && (rowIdx < n_subjects)){
	    atomicAdd(&d_score[pIdx], shared_score[0]);
	}
}


static __global__ void calculate_Ts(float * d_a, const float * d_b, float * d_Ts,
                       size_t n_permutations, size_t n_subjects, size_t voxel, bool * d_boolean_score) {
	if(d_boolean_score[voxel] == false)
	  return;
	const size_t pIdx = threadIdx.x + blockDim.x * blockIdx.x;
	if ((pIdx < (n_permutations + 1))) {
	    float l_a_1 = d_a[pIdx];
	    float l_a_2 = d_a[pIdx + (n_permutations + 1)];
/*
	   if(l_a_1 < 0.f)
	     l_a_1 = 0.f;

	    if(l_a_2 < 0.f)
	      l_a_2 = 0.f;*/


	    float row_0 = l_a_1*d_b[0] + l_a_2*d_b[1];


	    row_0 = row_0*row_0;
	    //row_0 = row_0*row_0;

	    float row_1 = l_a_1*d_b[2] + l_a_2*d_b[3];

	    row_1 = row_1*row_1;


	    d_Ts[pIdx] = 0.5f*(row_1 + row_0);
	}
}



static __global__ void sum_Ts(const float * d_Ts, const float * d_score, float * d_pval, size_t n_permutations, size_t voxel, bool * d_boolean_score){

  if(d_boolean_score[voxel] == false)
    return;

  const size_t rowIdx = threadIdx.x + blockIdx.x*blockDim.x;

  const size_t tIdx = threadIdx.x;

  extern __shared__ float shared_pval[];

  if (rowIdx < (n_permutations + 1)) {
      const float comp_val = d_Ts[0];
      if((d_Ts[rowIdx] >= comp_val  && d_score[rowIdx] >= 0.f)  || rowIdx == 0){
    	  shared_pval[tIdx] = 1.f;
      }else{
    	  shared_pval[tIdx] = 0.f;
      }
  }else{
      shared_pval[tIdx] = 0.f;

  }

  __syncthreads();

  for(size_t stride = blockDim.x/2 ; stride > 0 ; stride >>= 1){

      if(tIdx < stride  && ((rowIdx + stride) < (n_permutations + 1))){
    	  shared_pval[tIdx] += shared_pval[tIdx + stride];
      }

      __syncthreads();
  }

  if((threadIdx.x == 0) && (rowIdx < (n_permutations + 1))){
      atomicAdd(&d_pval[voxel], shared_pval[0]);
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
      d_pval[vIdx] = 1.f;
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


int compute_pvals(const float * d_sy, const float *  d_res,  const float * d_hat, const float * d_Sigma_A, const float * d_Sigma_E,
        float * d_pvals, const unsigned int * d_pmatrix,
        bool * d_boolean_score, aux_variables aux_vars, pval_variables pval_vars, bool covariate,
        size_t n_subjects, size_t n_voxels, size_t n_permutations, cudaStream_t stream, cublasHandle_t handle) {




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


	int blockSize_n_permutations = 32;

	if(n_permutations + 1 >= 64){
		blockSize_n_permutations = 64;
	}
	if(n_permutations + 1  >= 128){
		blockSize_n_permutations = 128;
	}
	if(n_permutations + 1  >= 256){
		blockSize_n_permutations = 256;
	}
	if(n_permutations + 1  >= 512){
		blockSize_n_permutations = 512;
	}
	if(n_permutations + 1  >= 1024){
		if((n_permutations + 1) % 1024 <= (n_permutations + 1) % 512){
			blockSize_n_permutations = 1024;
		}else{
			blockSize_n_permutations = 512;
		}

	}

	dim3 blockSize_1(blockSize_n_subjects,1, 1);
	dim3 gridSize_1(ceil(float(n_subjects)/float(blockSize_n_subjects)),n_permutations + 1, 1);

	dim3 blockSize_set(blockSize_n_subjects, 1, 1);

	dim3 gridSize_set(ceil(float(n_subjects)/float(blockSize_n_subjects)),ceil(float(n_permutations + 1)), 1);


	dim3 blockSize_2(blockSize_n_permutations, 1, 1);
	dim3 gridSize_2(ceil(float(n_permutations + 1)/float(blockSize_n_permutations)), 1, 1);



	 float alpha = 1.f;
	 float beta = 0.f;

	for(size_t voxel = 0 ; voxel < n_voxels ; voxel++){
		set_syP<<<gridSize_set, blockSize_set, 0, stream>>>(pval_vars.syP, (const float *)d_res, (const float *)d_sy, (const unsigned int *)d_pmatrix, voxel, n_permutations, n_subjects , n_voxels, d_boolean_score);
		gpuErrchk(cudaPeekAtLastError());
		if(covariate){
		    cublasErrchk(cublasSgemm_v2(handle, CUBLAS_OP_N, CUBLAS_OP_N, n_subjects,n_permutations + 1, n_subjects, &alpha,
		    		d_hat, n_subjects, pval_vars.syP, n_subjects, &beta, pval_vars.d_F, n_subjects));
		    set_F<<<gridSize_set, blockSize_set, 0, stream>>>(pval_vars.d_F,(const float*)pval_vars.d_F, d_pmatrix, n_subjects, n_permutations, voxel, d_boolean_score);
		    gpuErrchk(cudaPeekAtLastError());
		}else{
		    set_F<<<gridSize_set, blockSize_set, 0, stream>>>(pval_vars.d_F,(const float*)pval_vars.syP, d_pmatrix, n_subjects, n_permutations, voxel, d_boolean_score);
		    gpuErrchk(cudaPeekAtLastError());
		}



		gpuErrchk(cudaMemsetAsync(pval_vars.d_sigmaP, 0, sizeof(float) * (n_permutations + 1), stream));

		calculate_sigmaP<<<gridSize_1, blockSize_1, blockSize_n_subjects*sizeof(float), stream>>>((const float*)pval_vars.d_F, pval_vars.d_sigmaP, n_subjects, n_permutations, voxel, d_boolean_score);
		gpuErrchk(cudaPeekAtLastError());

		gpuErrchk(cudaMemsetAsync(pval_vars.d_score, 0, sizeof(float) * (n_permutations + 1), stream));

		calculate_score<<<gridSize_1, blockSize_1, blockSize_n_subjects*sizeof(float), stream>>>((const float*)aux_vars.d_Z,(const float*) pval_vars.d_F, (const float *)pval_vars.d_sigmaP,
		    pval_vars.d_score, n_subjects, n_permutations, voxel, d_boolean_score);
		gpuErrchk(cudaPeekAtLastError());
		update_syP<<<gridSize_set, blockSize_set, 0,stream>>>(pval_vars.d_F,(const float*) pval_vars.d_sigmaP, n_subjects, n_permutations, voxel,d_boolean_score);
		gpuErrchk(cudaPeekAtLastError());


	    cublasErrchk(cublasSgemm_v2(handle, CUBLAS_OP_T, CUBLAS_OP_N, n_permutations + 1, 2, n_subjects, &alpha, pval_vars.d_F, n_subjects,
	    		aux_vars.d_Z, n_subjects, &beta, pval_vars.d_a, n_permutations + 1));



		calculate_Ts<<<gridSize_2, blockSize_2, 0, stream>>>(pval_vars.d_a,(const float*) aux_vars.d_ZTZI, pval_vars.d_Ts,
				n_permutations, n_subjects, voxel, d_boolean_score);
		gpuErrchk(cudaPeekAtLastError());

		sum_Ts<<<gridSize_2, blockSize_2, blockSize_n_permutations*sizeof(float), stream>>>((const float *)pval_vars.d_Ts, pval_vars.d_score, d_pvals,  n_permutations, voxel,d_boolean_score);
		gpuErrchk(cudaPeekAtLastError());

	}
	mean_Ts<<<ceil(float(n_voxels)/float(BLOCK_SIZE_3)), BLOCK_SIZE_3, 0, stream >>>(d_pvals, n_permutations, n_voxels, d_boolean_score);
	gpuErrchk(cudaPeekAtLastError());


	return 1;

}
