project( MSc_Math_Finance_2018Doc)
cmake_minimum_required(VERSION 2.8)
OPTION(BUILDDOCUMENTATION "Build documentation" OFF)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)

IF(BUILDDOCUMENTATION)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
 add_custom_target(doc ALL
 ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
 COMMENT "Generating API documentation with Doxygen" VERBATIM 
 )
 
 add_custom_target(docpdf ALL
 make
 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc/latex
 COMMENT "Generating pdf" VERBATIM 
 )
 
ENDIF(BUILDDOCUMENTATION)

endif(DOXYGEN_FOUND)


