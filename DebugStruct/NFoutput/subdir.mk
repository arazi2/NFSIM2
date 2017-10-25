################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFoutput/NFoutput.cpp 

OBJS += \
./NFoutput/NFoutput.o 

CPP_DEPS += \
./NFoutput/NFoutput.d 


# Each subdirectory must supply rules for building sources it contributes
NFoutput/%.o: ../NFoutput/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


