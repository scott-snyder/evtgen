
########################################################################
# Copyright 1998-2020 CERN for the benefit of the EvtGen authors       #
#                                                                      #
# This file is part of EvtGen.                                         #
#                                                                      #
# EvtGen is free software: you can redistribute it and/or modify       #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# EvtGen is distributed in the hope that it will be useful,            #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with EvtGen.  If not, see <https://www.gnu.org/licenses/>.     #
########################################################################

# - Try to find Sherpa
# Defines:
#
#  SHERPA_FOUND
#  SHERPA_VERSION
#  SHERPA_INCLUDE_DIR
#  SHERPA_INCLUDE_DIRS (not cached)
#  SHERPA_LIBRARY
#  SHERPA_LIBRARIES (not cached)
#  SHERPA_LIBRARY_DIRS (not cached)

find_path(SHERPA_INCLUDE_DIR SHERPA/Main/Sherpa.H
          HINTS $ENV{SHERPA_ROOT_DIR}/include ${SHERPA_ROOT_DIR}/include
          PATH_SUFFIXES include SHERPA-MC)

mark_as_advanced(SHERPA_INCLUDE_DIR)

# Enforce a minimal list if none is explicitly requested
if(NOT Sherpa_FIND_COMPONENTS)
    set(Sherpa_FIND_COMPONENTS SherpaMain ModelMain METoolsMain METoolsLoops METoolsSpinCorrelations PDF PhotonsMEs PhotonsMain PhotonsPhaseSpace PhotonsTools Remnant ToolsMath ToolsOrg ToolsPhys)
endif()

foreach(component ${Sherpa_FIND_COMPONENTS})
  find_library(SHERPA_${component}_LIBRARY NAMES ${component} ${_${component}_names}
               HINTS $ENV{SHERPA_ROOT_DIR}/lib ${SHERPA_ROOT_DIR}/lib $ENV{SHERPA_ROOT_DIR}/lib64 ${SHERPA_ROOT_DIR}/lib64
               PATH_SUFFIXES SHERPA-MC)
  if (SHERPA_${component}_LIBRARY)
    set(SHERPA_${component}_FOUND 1)
    list(APPEND SHERPA_LIBRARIES ${SHERPA_${component}_LIBRARY})

    get_filename_component(libdir ${SHERPA_${component}_LIBRARY} PATH)
    list(APPEND SHERPA_LIBRARY_DIRS ${libdir})
  else()
    set(SHERPA_${component}_FOUND 0)
  endif()
  mark_as_advanced(SHERPA_${component}_LIBRARY)
endforeach()

if(SHERPA_LIBRARY_DIRS)
  list(REMOVE_DUPLICATES SHERPA_LIBRARY_DIRS)
endif()

find_program(SHERPA_CONFIG Sherpa-config
             HINTS $ENV{SHERPA_ROOT_DIR}/bin ${SHERPA_ROOT_DIR}/bin)

if(SHERPA_CONFIG)
    execute_process(COMMAND ${SHERPA_CONFIG} --version
                    OUTPUT_VARIABLE SHERPA_VERSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(SHERPA_VERSION ${SHERPA_VERSION} CACHE STRING "Detected version of Sherpa")
    mark_as_advanced(SHERPA_VERSION)
endif()

# handle the QUIETLY and REQUIRED arguments and set SHERPA_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Sherpa DEFAULT_MSG SHERPA_VERSION SHERPA_INCLUDE_DIR SHERPA_LIBRARIES )

set(SHERPA_INCLUDE_DIRS ${SHERPA_INCLUDE_DIR})

mark_as_advanced(SHERPA_FOUND)
