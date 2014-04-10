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
CMAKE_MINIMUM_REQUIRED (VERSION 2.8)
#message(${SIMS_SOURCE_DIR})
### DEPENDENCIES ########################################################################################

# essential dependencies
INCLUDE (${SIMS_SOURCE_DIR}/FindMECHSYS.cmake)
#message("idiota " ${MECHSYS_SOURCE_DIR})
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
