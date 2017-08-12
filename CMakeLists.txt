set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPolyFEMApplication ****")

######### add include directory here
include_directories( ${CMAKE_SOURCE_DIR}/external_libraries )
include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications )
#######################################################################

######### add additional cmake module path
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/applications/distributed_builders_application/custom_cmake_modules")
#######################################################################

######### add additional libraries
# IF(DEFINED PETSC_DIR)
#    SET(PETSC_EXECUTABLE_RUNS "YES") #skip the test
#    FIND_PACKAGE(PETSc REQUIRED)
#    IF(NOT ${PETSC_FOUND})
#        message(ERROR "PETSc is required for PartiallySaturatedSoilApplication")
#    ENDIF()
#    INCLUDE_DIRECTORIES(${PETSC_INCLUDES})
# ENDIF()


## generate variables with the sources
set( KRATOS_POLYFEM_APPLICATION_SOURCES     
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_quadrature/polygon_integration_rules.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/polyfem_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/polyfem_python_application.cpp
)


###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosPolyFEMApplication SHARED ${KRATOS_POLYFEM_APPLICATION_SOURCES})

target_link_libraries(KratosPolyFEMApplication KratosCore)
# target_link_libraries(KratosPolyFEMApplication ${PETSC_LIBRARIES})
# target_link_libraries(KratosPolyFEMApplication ${MPI_CXX_LIBRARY} ${MPI_C_LIBRARY})
set_target_properties(KratosPolyFEMApplication PROPERTIES PREFIX "")
install(TARGETS KratosPolyFEMApplication DESTINATION libs )

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosPolyFEMApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphysics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PolyFEMApplication.py" DESTINATION KratosMultiphysics )
