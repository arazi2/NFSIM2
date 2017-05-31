################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFcore/complex.cpp \
../NFcore/complexList.cpp \
../NFcore/molecule.cpp \
../NFcore/moleculeType.cpp \
../NFcore/observable.cpp \
../NFcore/reactionClass.cpp \
../NFcore/system.cpp \
../NFcore/templateMolecule.cpp 

OBJS += \
./NFcore/complex.o \
./NFcore/complexList.o \
./NFcore/molecule.o \
./NFcore/moleculeType.o \
./NFcore/observable.o \
./NFcore/reactionClass.o \
./NFcore/system.o \
./NFcore/templateMolecule.o 

CPP_DEPS += \
./NFcore/complex.d \
./NFcore/complexList.d \
./NFcore/molecule.d \
./NFcore/moleculeType.d \
./NFcore/observable.d \
./NFcore/reactionClass.d \
./NFcore/system.d \
./NFcore/templateMolecule.d 


# Each subdirectory must supply rules for building sources it contributes
NFcore/%.o: ../NFcore/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


