################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFinput/NFinput.cpp \
../NFinput/commandLineParser.cpp \
../NFinput/parseFuncXML.cpp \
../NFinput/parseSymRxns.cpp \
../NFinput/rnfRunner.cpp \
../NFinput/walk.cpp 

OBJS += \
./NFinput/NFinput.o \
./NFinput/commandLineParser.o \
./NFinput/parseFuncXML.o \
./NFinput/parseSymRxns.o \
./NFinput/rnfRunner.o \
./NFinput/walk.o 

CPP_DEPS += \
./NFinput/NFinput.d \
./NFinput/commandLineParser.d \
./NFinput/parseFuncXML.d \
./NFinput/parseSymRxns.d \
./NFinput/rnfRunner.d \
./NFinput/walk.d 


# Each subdirectory must supply rules for building sources it contributes
NFinput/%.o: ../NFinput/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


