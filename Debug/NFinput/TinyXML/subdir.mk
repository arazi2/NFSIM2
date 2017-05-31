################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFinput/TinyXML/tinystr.cpp \
../NFinput/TinyXML/tinyxml.cpp \
../NFinput/TinyXML/tinyxmlerror.cpp \
../NFinput/TinyXML/tinyxmlparser.cpp 

OBJS += \
./NFinput/TinyXML/tinystr.o \
./NFinput/TinyXML/tinyxml.o \
./NFinput/TinyXML/tinyxmlerror.o \
./NFinput/TinyXML/tinyxmlparser.o 

CPP_DEPS += \
./NFinput/TinyXML/tinystr.d \
./NFinput/TinyXML/tinyxml.d \
./NFinput/TinyXML/tinyxmlerror.d \
./NFinput/TinyXML/tinyxmlparser.d 


# Each subdirectory must supply rules for building sources it contributes
NFinput/TinyXML/%.o: ../NFinput/TinyXML/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


