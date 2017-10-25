################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFfunction/compositeFunction.cpp \
../NFfunction/funcParser.cpp \
../NFfunction/function.cpp \
../NFfunction/localFunction.cpp 

OBJS += \
./NFfunction/compositeFunction.o \
./NFfunction/funcParser.o \
./NFfunction/function.o \
./NFfunction/localFunction.o 

CPP_DEPS += \
./NFfunction/compositeFunction.d \
./NFfunction/funcParser.d \
./NFfunction/function.d \
./NFfunction/localFunction.d 


# Each subdirectory must supply rules for building sources it contributes
NFfunction/%.o: ../NFfunction/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


