################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFutil/conversion.cpp \
../NFutil/random.cpp \
../NFutil/stringOperations.cpp \
../NFutil/tools.cpp 

OBJS += \
./NFutil/conversion.o \
./NFutil/random.o \
./NFutil/stringOperations.o \
./NFutil/tools.o 

CPP_DEPS += \
./NFutil/conversion.d \
./NFutil/random.d \
./NFutil/stringOperations.d \
./NFutil/tools.d 


# Each subdirectory must supply rules for building sources it contributes
NFutil/%.o: ../NFutil/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


