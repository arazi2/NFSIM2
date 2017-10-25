################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFcore/moleculeLists/moleculeList.cpp 

OBJS += \
./NFcore/moleculeLists/moleculeList.o 

CPP_DEPS += \
./NFcore/moleculeLists/moleculeList.d 


# Each subdirectory must supply rules for building sources it contributes
NFcore/moleculeLists/%.o: ../NFcore/moleculeLists/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


