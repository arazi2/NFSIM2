################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include nauty24/subdir.mk
-include NFutil/MTrand/subdir.mk
-include NFutil/subdir.mk
-include NFtest/transcription/subdir.mk
-include NFtest/tlbr/subdir.mk
-include NFtest/simple_system/subdir.mk
-include NFtest/agentcell/cell/subdir.mk
-include NFtest/agentcell/subdir.mk
-include NFscheduler/subdir.mk
-include NFreactions/transformations/subdir.mk
-include NFreactions/reactions/subdir.mk
-include NFreactions/reactantLists/subdir.mk
-include NFreactions/mappings/subdir.mk
-include NFreactions/subdir.mk
-include NFoutput/subdir.mk
-include NFinput/TinyXML/subdir.mk
-include NFinput/subdir.mk
-include NFfunction/muParser/subdir.mk
-include NFfunction/subdir.mk
-include NFcore/reactionSelector/subdir.mk
-include NFcore/moleculeLists/subdir.mk
-include NFcore/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: src

# Tool invocations
src: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Cross G++ Linker'
	g++  -o "src" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(OBJS)$(C_UPPER_DEPS)$(CXX_DEPS)$(C_DEPS)$(CPP_DEPS) src
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
