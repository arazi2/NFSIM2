################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFcore/reactionSelector/directSelector.cpp \
../NFcore/reactionSelector/logClassSelector.cpp 

OBJS += \
./NFcore/reactionSelector/directSelector.o \
./NFcore/reactionSelector/logClassSelector.o 

CPP_DEPS += \
./NFcore/reactionSelector/directSelector.d \
./NFcore/reactionSelector/logClassSelector.d 


# Each subdirectory must supply rules for building sources it contributes
NFcore/reactionSelector/%.o: ../NFcore/reactionSelector/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


