################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFtest/simple_system/simple_system.cpp 

OBJS += \
./NFtest/simple_system/simple_system.o 

CPP_DEPS += \
./NFtest/simple_system/simple_system.d 


# Each subdirectory must supply rules for building sources it contributes
NFtest/simple_system/%.o: ../NFtest/simple_system/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


