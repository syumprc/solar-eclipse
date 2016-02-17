################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NiftiOrientConvert.cpp \
../src/RicVolume.cpp \
../src/RicVolumeSet.cpp 

OBJS += \
./src/NiftiOrientConvert.o \
./src/RicVolume.o \
./src/RicVolumeSet.o 

CPP_DEPS += \
./src/NiftiOrientConvert.d \
./src/RicVolume.d \
./src/RicVolumeSet.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I../../include -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


