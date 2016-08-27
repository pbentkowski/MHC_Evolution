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
CC=g++
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
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
	${OBJECTDIR}/Antigen.o \
	${OBJECTDIR}/DataHandler.o \
	${OBJECTDIR}/Environment.o \
	${OBJECTDIR}/Gene.o \
	${OBJECTDIR}/H2Pinteraction.o \
	${OBJECTDIR}/Host.o \
	${OBJECTDIR}/Pathogen.o \
	${OBJECTDIR}/RandomNumbs.o \
	${OBJECTDIR}/Tagging_system.o \
	${OBJECTDIR}/main.o


# C Compiler Flags
CFLAGS=-O3

# CC Compiler Flags
CCFLAGS=-m64 -O3 -std=c++1y
CXXFLAGS=-m64 -O3 -std=c++1y

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mhcevolution

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mhcevolution: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mhcevolution ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/Antigen.o: Antigen.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -w -I../build/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Antigen.o Antigen.cpp

${OBJECTDIR}/DataHandler.o: DataHandler.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -w -I../build/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/DataHandler.o DataHandler.cpp

${OBJECTDIR}/Environment.o: Environment.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -w -I../build/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Environment.o Environment.cpp

${OBJECTDIR}/Gene.o: Gene.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -w -I../build/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Gene.o Gene.cpp

${OBJECTDIR}/H2Pinteraction.o: H2Pinteraction.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -w -I../build/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/H2Pinteraction.o H2Pinteraction.cpp

${OBJECTDIR}/Host.o: Host.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -w -I../build/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Host.o Host.cpp

${OBJECTDIR}/Pathogen.o: Pathogen.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -w -I../build/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Pathogen.o Pathogen.cpp

${OBJECTDIR}/RandomNumbs.o: RandomNumbs.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -w -I../build/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/RandomNumbs.o RandomNumbs.cpp

${OBJECTDIR}/Tagging_system.o: Tagging_system.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -w -I../build/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Tagging_system.o Tagging_system.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -w -I../build/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mhcevolution

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
