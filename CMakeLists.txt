########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2009 Sergio Torres                                     #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# any later version.                                                   #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>  #
########################################################################

PROJECT                (SIMS)
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

# Flags
OPTION(A_MAKE_VERBOSE       "Show additional messages during compilation/link?" OFF)
OPTION(A_MAKE_ALL_WARNINGS  "Make with all warnings (-Wall)"                    OFF)
OPTION(A_MAKE_DEBUG_SYMBOLS "Make with debug symbols (-g)"                      OFF)
OPTION(A_MAKE_PROFILING     "Make with profiling information (-pg)"             OFF)
OPTION(A_MAKE_OPTIMIZED     "Make optimized (-O3)"                              ON)

ADD_DEFINITIONS (-fPIC -fmessage-length=0) # Each error message will appear on a single line; no line-wrapping will be done.
#INCLUDE_DIRECTORIES ($ENV{HOME}/pkg/boost_1_49_0)
### FLAGS ###############################################################################################

IF(A_MAKE_VERBOSE)
	SET(CMAKE_VERBOSE_MAKEFILE TRUE)
ENDIF(A_MAKE_VERBOSE)

IF(A_MAKE_ALL_WARNINGS)
	ADD_DEFINITIONS(-Wall)
ENDIF(A_MAKE_ALL_WARNINGS)

IF(A_MAKE_DEBUG_SYMBOLS)
	ADD_DEFINITIONS(-g)
ENDIF(A_MAKE_DEBUG_SYMBOLS)

IF(A_MAKE_PROFILING)
	SET(CMAKE_CXX_FLAGS -pg)
ENDIF(A_MAKE_PROFILING)

IF(A_MAKE_OPTIMIZED)
	ADD_DEFINITIONS(-Ofast -march=native)
ENDIF(A_MAKE_OPTIMIZED)

### DEPENDENCIES ########################################################################################

# essential dependencies
INCLUDE (${SIMS_SOURCE_DIR}/FindMECHSYS.cmake)
if(MECHSYS_FOUND)
    INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindDEPS.cmake)
    if(NOT DEPS_OK)
#        MESSAGE("Not all dependencies were satisfied")
    endif(NOT DEPS_OK)
    INCLUDE_DIRECTORIES (${MECHSYS_INCLUDE_DIRS})
else(MECHSYS_FOUND)
    MESSAGE("This project requires MechSys")
endif(MECHSYS_FOUND)

### SUB DIRECTORIES #####################################################################################

SUBDIRS(dem)
SUBDIRS(lbm)
