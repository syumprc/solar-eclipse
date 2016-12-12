/*
 * cudafphi.cuh
 *
 *  Created on: Aug 6, 2016
 *      Author: Brian Donohue
 */
#include <cublas_v2.h>
#include <stdio.h>
#ifndef CUDAFPHI_CUH_
#define CUDAFPHI_CUH_

static inline const char* cublasGetErrorString(cublasStatus_t status)
{
    switch(status)
    {
        case CUBLAS_STATUS_SUCCESS: return "CUBLAS_STATUS_SUCCESS";
        case CUBLAS_STATUS_NOT_INITIALIZED: return "CUBLAS_STATUS_NOT_INITIALIZED";
        case CUBLAS_STATUS_ALLOC_FAILED: return "CUBLAS_STATUS_ALLOC_FAILED";
        case CUBLAS_STATUS_INVALID_VALUE: return "CUBLAS_STATUS_INVALID_VALUE";
        case CUBLAS_STATUS_ARCH_MISMATCH: return "CUBLAS_STATUS_ARCH_MISMATCH";
        case CUBLAS_STATUS_MAPPING_ERROR: return "CUBLAS_STATUS_MAPPING_ERROR";
        case CUBLAS_STATUS_EXECUTION_FAILED: return "CUBLAS_STATUS_EXECUTION_FAILED";
        case CUBLAS_STATUS_INTERNAL_ERROR: return "CUBLAS_STATUS_INTERNAL_ERROR";
    }
    return "unknown error";
}


#define cublasErrchk(ans) { cublasAssert((ans), __FILE__, __LINE__); }
inline void cublasAssert(cublasStatus_t code, const char *file, int line, bool abort=true)
{

   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cublasGetErrorString(code), file, line);
      fprintf(stdout,"GPUassert: %s %s %d\n", cublasGetErrorString(code), file, line);

      if (abort) exit(code);
   }else{
	   fprintf(stderr,"GPUassert: %s %s %d\n", cublasGetErrorString(code), file, line);
   }
}

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{

   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      fprintf(stdout,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);

      if (abort) exit(code);
   }else{
	   fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
   }
}


#define BLOCK_SIZE_1 256
#define BLOCK_SIZE_2 512
#define BLOCK_SIZE_3 1024
#define BLOCK_SIZE_HAT 32
typedef struct {

   float * d_theta;
   float * d_weights;
   float * d_A, * d_B, * d_C, * d_D, * d_E;
   float * d_score;
   float * d_Sigma_A;
   float * d_Sigma_E;
   float * d_Sigma_P;

}compute_h2_variables;




typedef struct {
  float * mean_or_sigma;
  float * d_Y;
}compute_F_variables;


typedef struct {
  const float * d_Z;
  float * d_ZTZI;
}aux_variables;




typedef struct{
  float * syP;
  float * d_F;
  float * d_score;
  float * d_Ts;
  float * d_a;
  float * d_sigmaP;
}pval_variables;



typedef struct{

	float * hat;
	const unsigned int * pmatrix;
	const float * evectors;
	const float * Z;

	aux_variables aux_vars;

	bool covariates;
	bool get_pval;
	bool get_connectivity;

	int device_id;

	float * h2;
	float * indicator;
	float * pvals;



}cuda_fphi_variables_per_device;

typedef struct{

	cuda_fphi_variables_per_device shared_device_vars;
	int stream_number;

	int cov_number;
	size_t n_voxels;
	size_t n_subjects;
	size_t n_permutations;
	size_t n_covariates;
	float * h_y;

	float * h_cov;

}cuda_fphi_variables_per_stream;


int compute_h2(float * d_F, float * d_h2,  float * d_indicator, bool * d_boolean_score,
               compute_h2_variables vars, aux_variables aux_vars,
               size_t n_subjects, size_t n_voxels, cudaStream_t stream);


int compute_F(const float * d_hat, float* d_sy, const float *d_evectors, float * d_F,
              compute_F_variables vars, bool covariates, size_t n_subjects,
               size_t n_voxels, cudaStream_t stream, cublasHandle_t handle);

int compute_pvals(const float * d_sy, const float *  d_res,  const float * d_hat, const float * d_Sigma_A,
                  const float * d_sigma_E, float * d_pvals, const unsigned int * d_pmatrix,
                  bool * h_boolean_score, aux_variables aux_vars, pval_variables pval_vars,  bool covariate,
                  size_t n_subjects, size_t n_voxels, size_t n_permutations, cudaStream_t stream, cublasHandle_t handle);


void * run_cudafphi_pthread(void * cudafphi_args);
size_t run_allocation_test(cuda_fphi_variables_per_stream cudafPHI_args);
void * run_cudafphi_connect_pthread(void * cudafphi_args);

void whitened_Y(float * d_sy, float * mean, const float * d_evectors, float * d_y, size_t n_subjects, size_t n_voxels,cublasHandle_t handle, cudaStream_t stream);

#endif /* CUDAFPHI_H_ */
