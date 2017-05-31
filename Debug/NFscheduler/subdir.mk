################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFscheduler/NFstream.cpp \
../NFscheduler/Scheduler.cpp 

OBJS += \
./NFscheduler/NFstream.o \
./NFscheduler/Scheduler.o 

CPP_DEPS += \
./NFscheduler/NFstream.d \
./NFscheduler/Scheduler.d 


# Each subdirectory must supply rules for building sources it contributes
NFscheduler/%.o: ../NFscheduler/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


