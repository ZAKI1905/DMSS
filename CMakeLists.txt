cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(DMSS VERSION 1.0 DESCRIPTION "DMSS: Dark Matter Solar Signals")

string(TIMESTAMP DMSS_RELEASE_DATE "%m, %d, %Y")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

set(DEP_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/lib)

cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)

if(${_host_name} STREQUAL "Mohammadrezas-MacBook-Pro.local")
    set(ROOT_LIB_DIR /opt/local/lib/root6)
    set(ROOT_INC_DIR /opt/local/libexec/root6/include/root)
    set(MINUIT_INC_DIR /opt/local/include/root6/Minuit2)
    set(MINUIT_LIB_DIR /opt/local/lib/root6)
    set(OS_NAME darwin)

else()
    message("Error: Hostname unknown!")
    message("Error: Setup the include and library paths in 'CMakeLists.txt'.")
    message("Your host name: " ${_host_name})
endif()

set(DEP_INC_DIR dependencies/include)

find_package(GSL REQUIRED)

set(ROOT_LIBS
    ${ROOT_LIB_DIR}/libCore.so
    ${ROOT_LIB_DIR}/libGpad.so
    ${ROOT_LIB_DIR}/libGraf.so
    ${ROOT_LIB_DIR}/libHist.so
    ${ROOT_LIB_DIR}/libMathCore.so
    ${ROOT_LIB_DIR}/libMathMore.so
    ${ROOT_LIB_DIR}/libGenVector.so
    ${ROOT_LIB_DIR}/libPhysics.so
    ${ROOT_LIB_DIR}/libRooFit.so
    ${ROOT_LIB_DIR}/libRooFitCore.so
    ${ROOT_LIB_DIR}/libRooStats.so
    ${ROOT_LIB_DIR}/libRIO.so
)

set(MINUIT_LIB ${MINUIT_LIB_DIR}/libMinuit2.so)
set(ZAKI_LIB ${DEP_LIB_DIR}/Zaki/${OS_NAME}/libZaki.a)
set(CONFIND_LIB ${DEP_LIB_DIR}/Confind/${OS_NAME}/libConfind.a)

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
                            ${MINUIT_INC_DIR}
                            ${GSL_INCLUDE_DIRS}
                            ${ROOT_INC_DIR}
                            ${DEP_INC_DIR}
                            )

target_link_libraries(DMSS PUBLIC 
                        ${ROOT_LIBS} 
                        ${ZAKI_LIB} 
                        ${CONFIND_LIB} 
                        ${MINUIT_LIB}
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
                                ${MINUIT_INC_DIR}
                                ${GSL_INC_DIR}
                                ${ROOT_INC_DIR}
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