# ...........................................................
#                   EDIT THIS SECTION ONLY
# ...........................................................
# Add the path to "ROOTConfig.cmake" below :

set(root_cmake_path "/opt/local/libexec/root6/share/root/cmake")

# ...........................................................
#          !    DON'T EDIT BELOW THIS LINE    !
# ...........................................................
cmake_minimum_required(VERSION 3.10)

project(DMSS VERSION 1.0 DESCRIPTION "DMSS: Dark Matter Solar Signals")


string(TIMESTAMP DMSS_RELEASE_DATE "%m, %d, %Y")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

set(DEP_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/lib)
set(DEP_INC_DIR dependencies/include)

set(CMAKE_PREFIX_PATH ${root_cmake_path})

find_package(GSL REQUIRED)
find_package(ROOT COMPONENTS Minuit2)

set(ROOT_LIBS
    ROOT::Core
    ROOT::Gpad
    ROOT::Graf
    ROOT::Hist
    ROOT::MathCore
    ROOT::MathMore
    ROOT::GenVector
    ROOT::Physics
    ROOT::RooFit
    ROOT::RooFitCore
    ROOT::RooStats
    ROOT::RIO
)

set(ZAKI_LIB ${DEP_LIB_DIR}/Zaki/${CMAKE_SYSTEM_NAME}/libZaki.a)
set(CONFIND_LIB ${DEP_LIB_DIR}/Confind/${CMAKE_SYSTEM_NAME}/libConfind.a)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file(DMSSConfig.h.in ${CMAKE_CURRENT_SOURCE_DIR}/include/DMSS/DMSSConfig.h)

find_package(OpenMP REQUIRED)
IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

add_subdirectory(include)
add_subdirectory(src)

add_library(DMSS STATIC ${DMSS_SRC_Files})
set_target_properties(DMSS PROPERTIES VERSION ${PROJECT_VERSION})

target_include_directories(DMSS 
                            PUBLIC
                            ${CMAKE_CURRENT_SOURCE_DIR}/include
                            ${GSL_INCLUDE_DIRS}
                            ${ROOT_INCLUDE_DIRS}
                            ${DEP_INC_DIR}
                            )

target_link_libraries(DMSS PUBLIC 
                        ${ROOT_LIBS} 
                        ${ZAKI_LIB} 
                        ${CONFIND_LIB} 
                        ROOT::Minuit2
                        GSL::gsl
                        GSL::gslcblas
                        OpenMP::OpenMP_CXX) 

                            
install(TARGETS DMSS DESTINATION lib)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/include"
        DESTINATION dependencies)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/lib" 
        DESTINATION dependencies)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/License" DESTINATION Docs)

# Option for generating the doxygen docs
option(GEN_DOCS "Generate documents using doxygen." OFF)

if(GEN_DOCS)
    execute_process (
        COMMAND
        doxygen ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

add_subdirectory(main)

foreach(src ${DMSS_main_sources})
    add_executable(${src} "main/${src}.cpp")  
    target_link_libraries(${src} PUBLIC DMSS)
    target_include_directories(${src} PUBLIC  
                                ${CMAKE_CURRENT_SOURCE_DIR}/include
                                ${GSL_INCLUDE_DIRS}
                                ${ROOT_INCLUDE_DIRS}
                                ${DEP_INC_DIR}
                                )
    install(FILES "main/${src}.cpp" DESTINATION main)
    install(TARGETS ${src} DESTINATION main/${src})
endforeach()

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${DMSS_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${DMSS_VERSION_MINOR}")
include(CPack)