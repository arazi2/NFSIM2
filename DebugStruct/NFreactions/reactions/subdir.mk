################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFreactions/reactions/DORreaction.cpp \
../NFreactions/reactions/reaction.cpp \
../NFreactions/reactions/RHSreaction.cpp 

OBJS += \
./NFreactions/reactions/DORreaction.o \
./NFreactions/reactions/reaction.o \
./NFreactions/reactions/RHSreaction.o 


CPP_DEPS += \
./NFreactions/reactions/DORreaction.d \
./NFreactions/reactions/reaction.d \
./NFreactions/reactions/RHSreaction.d 


# Each subdirectory must supply rules for building sources it contributes
NFreactions/reactions/%.o: ../NFreactions/reactions/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


