
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


# the test executables require ROOT for histogramming etc.
find_package(ROOT 6.18 QUIET COMPONENTS Hist)
find_package(nlohmann_json 3.6.1 QUIET)
find_package(Threads QUIET)
find_package(TBB QUIET)
if (ROOT_FOUND AND nlohmann_json_FOUND AND Threads_FOUND AND TBB_FOUND)
    # build each of the executables
    foreach( test_exe evt_dalitz evtgenlhc_test1 example1 exampleWriteHepMC testDecayModel )
        add_executable(${test_exe} ${test_exe}.cc)
        target_link_libraries(${test_exe} PRIVATE EvtGen)
        if (EVTGEN_HEPMC3)
            target_compile_definitions(${test_exe} PRIVATE EVTGEN_HEPMC3)
        endif()
        if( EVTGEN_PYTHIA OR EVTGEN_PHOTOS OR EVTGEN_TAUOLA OR EVTGEN_SHERPA )
            target_compile_definitions(${test_exe} PRIVATE EVTGEN_EXTERNAL)
            target_link_libraries(${test_exe} PRIVATE EvtGenExternal)
        endif()
        target_link_libraries(${test_exe} PRIVATE ROOT::Core ROOT::RIO ROOT::Hist ROOT::Tree)
    endforeach()

    # testDecayModel needs nlohmann_json and Threads as well
    target_link_libraries(testDecayModel PRIVATE nlohmann_json::nlohmann_json)
    target_link_libraries(testDecayModel PRIVATE Threads::Threads)
    target_link_libraries(testDecayModel PRIVATE TBB::tbb)

    # install the executables
    install(TARGETS evt_dalitz evtgenlhc_test1 example1 exampleWriteHepMC testDecayModel
        DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/test
        )

    # install the decay files, macros, scripts, etc.
    install(DIRECTORY exampleFiles DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/test)
    install(DIRECTORY jsonFiles DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/test)
    install(FILES do_tests runTests.py makeSrcDepsJson makeModelsJson
        DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/test
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
else()
    message(WARNING "Could not find ROOT and/or nlohmann_json, cannot build test executables")
endif()
