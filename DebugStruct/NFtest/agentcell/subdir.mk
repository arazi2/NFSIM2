################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFtest/agentcell/agentcell.cpp 

OBJS += \
./NFtest/agentcell/agentcell.o 

CPP_DEPS += \
./NFtest/agentcell/agentcell.d 


# Each subdirectory must supply rules for building sources it contributes
NFtest/agentcell/%.o: ../NFtest/agentcell/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


