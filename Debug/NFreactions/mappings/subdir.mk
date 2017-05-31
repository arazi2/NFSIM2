################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFreactions/mappings/mapping.cpp \
../NFreactions/mappings/mappingGenerator.cpp \
../NFreactions/mappings/mappingSet.cpp 

OBJS += \
./NFreactions/mappings/mapping.o \
./NFreactions/mappings/mappingGenerator.o \
./NFreactions/mappings/mappingSet.o 

CPP_DEPS += \
./NFreactions/mappings/mapping.d \
./NFreactions/mappings/mappingGenerator.d \
./NFreactions/mappings/mappingSet.d 


# Each subdirectory must supply rules for building sources it contributes
NFreactions/mappings/%.o: ../NFreactions/mappings/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


