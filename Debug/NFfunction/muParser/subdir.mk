################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../NFfunction/muParser/muParser.cpp \
../NFfunction/muParser/muParserBase.cpp \
../NFfunction/muParser/muParserBytecode.cpp \
../NFfunction/muParser/muParserCallback.cpp \
../NFfunction/muParser/muParserComplex.cpp \
../NFfunction/muParser/muParserError.cpp \
../NFfunction/muParser/muParserInt.cpp \
../NFfunction/muParser/muParserTokenReader.cpp 

OBJS += \
./NFfunction/muParser/muParser.o \
./NFfunction/muParser/muParserBase.o \
./NFfunction/muParser/muParserBytecode.o \
./NFfunction/muParser/muParserCallback.o \
./NFfunction/muParser/muParserComplex.o \
./NFfunction/muParser/muParserError.o \
./NFfunction/muParser/muParserInt.o \
./NFfunction/muParser/muParserTokenReader.o 

CPP_DEPS += \
./NFfunction/muParser/muParser.d \
./NFfunction/muParser/muParserBase.d \
./NFfunction/muParser/muParserBytecode.d \
./NFfunction/muParser/muParserCallback.d \
./NFfunction/muParser/muParserComplex.d \
./NFfunction/muParser/muParserError.d \
./NFfunction/muParser/muParserInt.d \
./NFfunction/muParser/muParserTokenReader.d 


# Each subdirectory must supply rules for building sources it contributes
NFfunction/muParser/%.o: ../NFfunction/muParser/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


