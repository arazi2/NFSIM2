################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFreactions/transformations/moleculeCreator.cpp \
../NFreactions/transformations/speciesCreator.cpp \
../NFreactions/transformations/transformation.cpp \
../NFreactions/transformations/transformationSet.cpp 

OBJS += \
./NFreactions/transformations/moleculeCreator.o \
./NFreactions/transformations/speciesCreator.o \
./NFreactions/transformations/transformation.o \
./NFreactions/transformations/transformationSet.o 

CPP_DEPS += \
./NFreactions/transformations/moleculeCreator.d \
./NFreactions/transformations/speciesCreator.d \
./NFreactions/transformations/transformation.d \
./NFreactions/transformations/transformationSet.d 


# Each subdirectory must supply rules for building sources it contributes
NFreactions/transformations/%.o: ../NFreactions/transformations/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


