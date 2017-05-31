################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFreactions/reactantLists/reactantList.cpp \
../NFreactions/reactantLists/reactantTree.cpp 

OBJS += \
./NFreactions/reactantLists/reactantList.o \
./NFreactions/reactantLists/reactantTree.o 

CPP_DEPS += \
./NFreactions/reactantLists/reactantList.d \
./NFreactions/reactantLists/reactantTree.d 


# Each subdirectory must supply rules for building sources it contributes
NFreactions/reactantLists/%.o: ../NFreactions/reactantLists/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


