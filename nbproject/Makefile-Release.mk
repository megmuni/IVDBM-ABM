#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/examples/ex_diffuse_temperature.o \
	${OBJECTDIR}/examples/test_not_main.o \
	${OBJECTDIR}/findCudaTests/July3CUDA.o \
	${OBJECTDIR}/findCudaTests/July3main.o \
	${OBJECTDIR}/findCudaTests/conflict-main.o \
	${OBJECTDIR}/findCudaTests/conflict.o \
	${OBJECTDIR}/findCudaTests/conflict.o \
	${OBJECTDIR}/findCudaTests/main-ptx.o \
	${OBJECTDIR}/findCudaTests/main.o \
	${OBJECTDIR}/findCudaTests/main_for_lib.o \
	${OBJECTDIR}/findCudaTests/myTest.o \
	${OBJECTDIR}/findCudaTests/path\ with\ spaces/conflict.o \
	${OBJECTDIR}/findCudaTests/path\ with\ spaces/conflict.o \
	${OBJECTDIR}/findCudaTests/path\ with\ spaces/no-conflict.o \
	${OBJECTDIR}/findCudaTests/path\ with\ spaces/no-conflict.o \
	${OBJECTDIR}/findCudaTests/test_bin.o \
	${OBJECTDIR}/findCudaTests/test_lib.o \
	${OBJECTDIR}/findCudaTests/test_ptx.o \
	${OBJECTDIR}/findCudaTests/test_with_spaces.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/src/Agent/Agent.o \
	${OBJECTDIR}/src/Agent/Usr_Agents/Collagen.o \
	${OBJECTDIR}/src/Agent/Usr_Agents/Chondrocyte.o \
	${OBJECTDIR}/src/Agent/Usr_Agents/Hyaluronan.o \
	${OBJECTDIR}/src/Agent/Usr_Agents/Macrophage.o \
	${OBJECTDIR}/src/Agent/Usr_Agents/Neutrophil.o \
	${OBJECTDIR}/src/Agent/Usr_Agents/Platelet.o \
	${OBJECTDIR}/src/FieldVariable/FieldVariable.o \
	${OBJECTDIR}/src/FieldVariable/Usr_FieldVariables/WHChemical.o \
	${OBJECTDIR}/src/Patch/Patch.o \
	${OBJECTDIR}/src/World/Lattice/Lattice.o \
	${OBJECTDIR}/src/World/Usr_World/woundHealingWorld.o \
	${OBJECTDIR}/src/World/World.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/cpu_abm_apr14

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/cpu_abm_apr14: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/cpu_abm_apr14 ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/examples/ex_diffuse_temperature.o: examples/ex_diffuse_temperature.cpp 
	${MKDIR} -p ${OBJECTDIR}/examples
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/examples/ex_diffuse_temperature.o examples/ex_diffuse_temperature.cpp

${OBJECTDIR}/examples/test_not_main.o: examples/test_not_main.cpp 
	${MKDIR} -p ${OBJECTDIR}/examples
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/examples/test_not_main.o examples/test_not_main.cpp

${OBJECTDIR}/findCudaTests/July3CUDA.o: findCudaTests/July3CUDA.cu 
	${MKDIR} -p ${OBJECTDIR}/findCudaTests
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/findCudaTests/July3CUDA.o findCudaTests/July3CUDA.cu

${OBJECTDIR}/findCudaTests/July3main.o: findCudaTests/July3main.cu 
	${MKDIR} -p ${OBJECTDIR}/findCudaTests
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/findCudaTests/July3main.o findCudaTests/July3main.cu

${OBJECTDIR}/findCudaTests/conflict-main.o: findCudaTests/conflict-main.cpp 
	${MKDIR} -p ${OBJECTDIR}/findCudaTests
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/findCudaTests/conflict-main.o findCudaTests/conflict-main.cpp

${OBJECTDIR}/findCudaTests/conflict.o: findCudaTests/conflict.cpp 
	${MKDIR} -p ${OBJECTDIR}/findCudaTests
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/findCudaTests/conflict.o findCudaTests/conflict.cpp

${OBJECTDIR}/findCudaTests/conflict.o: findCudaTests/conflict.cu 
	${MKDIR} -p ${OBJECTDIR}/findCudaTests
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/findCudaTests/conflict.o findCudaTests/conflict.cu

${OBJECTDIR}/findCudaTests/main-ptx.o: findCudaTests/main-ptx.cpp 
	${MKDIR} -p ${OBJECTDIR}/findCudaTests
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/findCudaTests/main-ptx.o findCudaTests/main-ptx.cpp

${OBJECTDIR}/findCudaTests/main.o: findCudaTests/main.cc 
	${MKDIR} -p ${OBJECTDIR}/findCudaTests
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/findCudaTests/main.o findCudaTests/main.cc

${OBJECTDIR}/findCudaTests/main_for_lib.o: findCudaTests/main_for_lib.cc 
	${MKDIR} -p ${OBJECTDIR}/findCudaTests
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/findCudaTests/main_for_lib.o findCudaTests/main_for_lib.cc

${OBJECTDIR}/findCudaTests/myTest.o: findCudaTests/myTest.cu 
	${MKDIR} -p ${OBJECTDIR}/findCudaTests
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/findCudaTests/myTest.o findCudaTests/myTest.cu

.NO_PARALLEL:${OBJECTDIR}/findCudaTests/path\ with\ spaces/conflict.o
${OBJECTDIR}/findCudaTests/path\ with\ spaces/conflict.o: findCudaTests/path\ with\ spaces/conflict.cpp 
	${MKDIR} -p ${OBJECTDIR}/findCudaTests spaces
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/findCudaTests/path\ with\ spaces/conflict.o findCudaTests/path\ with\ spaces/conflict.cpp

.NO_PARALLEL:${OBJECTDIR}/findCudaTests/path\ with\ spaces/conflict.o
${OBJECTDIR}/findCudaTests/path\ with\ spaces/conflict.o: findCudaTests/path\ with\ spaces/conflict.cu 
	${MKDIR} -p ${OBJECTDIR}/findCudaTests spaces
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/findCudaTests/path\ with\ spaces/conflict.o findCudaTests/path\ with\ spaces/conflict.cu

.NO_PARALLEL:${OBJECTDIR}/findCudaTests/path\ with\ spaces/no-conflict.o
${OBJECTDIR}/findCudaTests/path\ with\ spaces/no-conflict.o: findCudaTests/path\ with\ spaces/no-conflict.cpp 
	${MKDIR} -p ${OBJECTDIR}/findCudaTests spaces
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/findCudaTests/path\ with\ spaces/no-conflict.o findCudaTests/path\ with\ spaces/no-conflict.cpp

.NO_PARALLEL:${OBJECTDIR}/findCudaTests/path\ with\ spaces/no-conflict.o
${OBJECTDIR}/findCudaTests/path\ with\ spaces/no-conflict.o: findCudaTests/path\ with\ spaces/no-conflict.cu 
	${MKDIR} -p ${OBJECTDIR}/findCudaTests spaces
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/findCudaTests/path\ with\ spaces/no-conflict.o findCudaTests/path\ with\ spaces/no-conflict.cu

${OBJECTDIR}/findCudaTests/test_bin.o: findCudaTests/test_bin.cu 
	${MKDIR} -p ${OBJECTDIR}/findCudaTests
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/findCudaTests/test_bin.o findCudaTests/test_bin.cu

${OBJECTDIR}/findCudaTests/test_lib.o: findCudaTests/test_lib.cu 
	${MKDIR} -p ${OBJECTDIR}/findCudaTests
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/findCudaTests/test_lib.o findCudaTests/test_lib.cu

${OBJECTDIR}/findCudaTests/test_ptx.o: findCudaTests/test_ptx.cu 
	${MKDIR} -p ${OBJECTDIR}/findCudaTests
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/findCudaTests/test_ptx.o findCudaTests/test_ptx.cu

${OBJECTDIR}/findCudaTests/test_with_spaces.o: findCudaTests/test_with_spaces.cpp 
	${MKDIR} -p ${OBJECTDIR}/findCudaTests
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/findCudaTests/test_with_spaces.o findCudaTests/test_with_spaces.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/src/Agent/Agent.o: src/Agent/Agent.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/Agent
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Agent/Agent.o src/Agent/Agent.cpp

${OBJECTDIR}/src/Agent/Usr_Agents/Collagen.o: src/Agent/Usr_Agents/Collagen.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/Agent/Usr_Agents
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Agent/Usr_Agents/Collagen.o src/Agent/Usr_Agents/Collagen.cpp

${OBJECTDIR}/src/Agent/Usr_Agents/Chondrocyte.o: src/Agent/Usr_Agents/Chondrocyte.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/Agent/Usr_Agents
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Agent/Usr_Agents/Chondrocyte.o src/Agent/Usr_Agents/Chondrocyte.cpp

${OBJECTDIR}/src/Agent/Usr_Agents/Hyaluronan.o: src/Agent/Usr_Agents/Hyaluronan.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/Agent/Usr_Agents
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Agent/Usr_Agents/Hyaluronan.o src/Agent/Usr_Agents/Hyaluronan.cpp

${OBJECTDIR}/src/Agent/Usr_Agents/Macrophage.o: src/Agent/Usr_Agents/Macrophage.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/Agent/Usr_Agents
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Agent/Usr_Agents/Macrophage.o src/Agent/Usr_Agents/Macrophage.cpp

${OBJECTDIR}/src/Agent/Usr_Agents/Neutrophil.o: src/Agent/Usr_Agents/Neutrophil.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/Agent/Usr_Agents
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Agent/Usr_Agents/Neutrophil.o src/Agent/Usr_Agents/Neutrophil.cpp

${OBJECTDIR}/src/Agent/Usr_Agents/Platelet.o: src/Agent/Usr_Agents/Platelet.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/Agent/Usr_Agents
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Agent/Usr_Agents/Platelet.o src/Agent/Usr_Agents/Platelet.cpp

${OBJECTDIR}/src/FieldVariable/FieldVariable.o: src/FieldVariable/FieldVariable.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/FieldVariable
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/FieldVariable/FieldVariable.o src/FieldVariable/FieldVariable.cpp

${OBJECTDIR}/src/FieldVariable/Usr_FieldVariables/WHChemical.o: src/FieldVariable/Usr_FieldVariables/WHChemical.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/FieldVariable/Usr_FieldVariables
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/FieldVariable/Usr_FieldVariables/WHChemical.o src/FieldVariable/Usr_FieldVariables/WHChemical.cpp

${OBJECTDIR}/src/Patch/Patch.o: src/Patch/Patch.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/Patch
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Patch/Patch.o src/Patch/Patch.cpp

${OBJECTDIR}/src/World/Lattice/Lattice.o: src/World/Lattice/Lattice.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/World/Lattice
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/World/Lattice/Lattice.o src/World/Lattice/Lattice.cpp

${OBJECTDIR}/src/World/Usr_World/woundHealingWorld.o: src/World/Usr_World/woundHealingWorld.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/World/Usr_World
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/World/Usr_World/woundHealingWorld.o src/World/Usr_World/woundHealingWorld.cpp

${OBJECTDIR}/src/World/World.o: src/World/World.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/World
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/World/World.o src/World/World.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/cpu_abm_apr14

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
