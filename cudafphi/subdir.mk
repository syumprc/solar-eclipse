################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../call_cudafPHI.cu \
../compute_F.cu \
../compute_h2.cu \
../compute_pval.cu \
../cudafPHI.cu \
../cudafphi_device_functions.cu 

CC_SRCS += \
../cuda_fphimain.cc 

CC_DEPS += \
./cuda_fphimain.d 

OBJS += \
./call_cudafPHI.o \
./compute_F.o \
./compute_h2.o \
./compute_pval.o \
./cuda_fphimain.o \
./cudafPHI.o \
./cudafphi_device_functions.o 

CU_DEPS += \
./call_cudafPHI.d \
./compute_F.d \
./compute_h2.d \
./compute_pval.d \
./cudafPHI.d \
./cudafphi_device_functions.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-7.5/bin/nvcc -G -g -O3 --use_fast_math -Xcompiler -fPIC -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_53,code=sm_53  -odir "." -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-7.5/bin/nvcc -G -g -O3 --use_fast_math -Xcompiler -fPIC --compile --relocatable-device-code=true -gencode arch=compute_35,code=compute_35 -gencode arch=compute_37,code=compute_37 -gencode arch=compute_50,code=compute_50 -gencode arch=compute_52,code=compute_52 -gencode arch=compute_53,code=compute_53 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_53,code=sm_53  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

%.o: ../%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-7.5/bin/nvcc -G -g -O3 --use_fast_math -Xcompiler -fPIC -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_53,code=sm_53  -odir "." -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-7.5/bin/nvcc -G -g -O3 --use_fast_math -Xcompiler -fPIC --compile  -x c++ -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


