################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../nauty24/nausparse.c \
../nauty24/nautil.c \
../nauty24/nauty.c 

OBJS += \
./nauty24/nausparse.o \
./nauty24/nautil.o \
./nauty24/nauty.o 

C_DEPS += \
./nauty24/nausparse.d \
./nauty24/nautil.d \
./nauty24/nauty.d 


# Each subdirectory must supply rules for building sources it contributes
nauty24/%.o: ../nauty24/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


