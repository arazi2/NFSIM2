################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFtest/transcription/transcription.cpp 

OBJS += \
./NFtest/transcription/transcription.o 

CPP_DEPS += \
./NFtest/transcription/transcription.d 


# Each subdirectory must supply rules for building sources it contributes
NFtest/transcription/%.o: ../NFtest/transcription/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


