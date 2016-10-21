/*
 * cudafPHI.cuh
 *
 *  Created on: Aug 6, 2016
 *      Author: brian
 */
#include <vector>


#ifndef CUDAFPHI_CUH_
#define CUDAFPHI_CUH_

#define BLOCK_SIZE_1 256
#define BLOCK_SIZE_2 1024
#define BLOCK_SIZE_3 512
#define BLOCK_SIZE_HAT 16
typedef struct {

   float * d_theta;
   float * d_weights;
   float * d_A, * d_B, * d_C, * d_D, * d_E;
   float * d_score;
   float * d_Sigma_P;

}compute_h2_variables;




typedef struct {
  float * mean;
  float * d_Y;
}compute_F_variables;

typedef struct{
  float * h_y;
  float * d_hat;
  float * d_sy;
  float * d_evectors;
  float * d_F;
  compute_F_variables compute_F_vars;
  bool covariates;
  size_t n_subjects;
  int current_stream_number;
  cudaStream_t stream;
}run_compute_F_variables;


typedef struct {
  float * d_Z;
  float * d_ZTZI;
}aux_variables;


typedef struct{
  float * d_F;
  float * d_h2;
  float * d_indicator;
  float * d_boolean_score;
  float * indicator;
  float * h2;
  compute_h2_variables compute_h2_vars;
  aux_variables aux_vars;
  size_t n_subjects;
  int stream_number;
  cudaStream_t stream;
}run_compute_h2_variables;



typedef struct{
  float * syP;
  float * d_F;
  float * d_score;
  float * d_Ts;
  float * d_a;
  float * d_sigmaP;
}pval_variables;


typedef struct{
  float * d_F;
  float * d_pvals;
  float * d_boolean_score;
  float * pvals;
  float * d_hat;
  compute_h2_variables compute_h2_vars;
  aux_variables aux_vars;
  pval_variables pval_vars;
  size_t n_subjects;
  size_t n_permutations;
  int stream_number;
  bool covariates;
  cudaStream_t stream;
}run_compute_pval_variables;

typedef struct{
  float * h_y;
  float * h2;
  float * indicator;
  float * pvals;
  float * d_sy;
  float * d_F;
  float * d_hat;
  float * d_evectors;
  float * d_h2;
  float * d_indicator;
  float * d_pvals;
  float * h_evectors;
  float * h_hat;
  float * h_Z;
  unsigned int * h_pmatrix;
  unsigned int * d_pmatrix;
  compute_F_variables compute_F_vars;
  aux_variables aux_vars;
  compute_h2_variables compute_h2_vars;
  pval_variables pval_vars;
  bool covariates;
  bool get_pval;
  bool * d_boolean_score;
  size_t n_voxels;
  size_t n_subjects;
  size_t n_permutations;
//  size_t voxels_used;
  int device;
  cudaStream_t stream;
  int stream_number;
}cudafPHI_variables;


int compute_h2(float * d_F, float * d_h2,  float * d_indicator, bool * d_boolean_score,
               compute_h2_variables vars, aux_variables aux_vars,
               size_t n_subjects, size_t n_voxels, cudaStream_t stream);


int compute_F(float * d_hat, float* d_sy, float *d_evectors, float * d_F,
              compute_F_variables vars, bool covariates, size_t n_subjects,
               size_t n_voxels, cudaStream_t stream);

int compute_pvals(float * d_sy, float * d_res, float * d_hat, float * d_sigma_E,
                  float * d_sigma_A, float * d_pvals, unsigned int * d_pmatrix,
                  bool * h_boolean_score, aux_variables aux_vars, pval_variables pval_vars, bool covariate,
                  size_t n_subjects, size_t n_voxels, size_t n_permutations, cudaStream_t stream);

void * run_cudafPHI_pthread(void * cudafPHI_args);
void  run_cudafPHI_loop(void * cudafPHI_args);
void * run_compute_F_pthread(void * run_compute_F_args);

void * run_compute_h2_pthread(void * run_compute_F_args);

void * run_compute_pval_pthread(void * run_compute_F_args);


void * run_allocation_test(void * cudafPHI_args);


#endif /* CUDAFPHI_H_ */
