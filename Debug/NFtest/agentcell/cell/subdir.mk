################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFtest/agentcell/cell/cell.cpp \
../NFtest/agentcell/cell/environment.cpp \
../NFtest/agentcell/cell/util.cpp 

OBJS += \
./NFtest/agentcell/cell/cell.o \
./NFtest/agentcell/cell/environment.o \
./NFtest/agentcell/cell/util.o 

CPP_DEPS += \
./NFtest/agentcell/cell/cell.d \
./NFtest/agentcell/cell/environment.d \
./NFtest/agentcell/cell/util.d 


# Each subdirectory must supply rules for building sources it contributes
NFtest/agentcell/cell/%.o: ../NFtest/agentcell/cell/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


