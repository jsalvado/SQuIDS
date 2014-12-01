i3_project(SQuIDS)

i3_add_library(SQuIDS
  private/SQuIDS/*.cpp
  USE_TOOLS gsl
  WITHOUT_I3_HEADERS)
set_source_files_properties(${PROJECT_SOURCE_DIR}/private/${PROJECT_NAME}/const.cpp  COMPILE_FLAGS "-std=c++11")
set_source_files_properties(${PROJECT_SOURCE_DIR}/private/${PROJECT_NAME}/SQUIDS.cpp COMPILE_FLAGS "-std=c++11")
set_source_files_properties(${PROJECT_SOURCE_DIR}/private/${PROJECT_NAME}/SUNalg.cpp COMPILE_FLAGS "-std=c++11 -Warray-bounds")

enable_testing()
file(GLOB tests ${PROJECT_SOURCE_DIR}/private/test/*.test.cpp)
foreach(t ${tests})
  get_filename_component(f ${t} NAME_WE)
  i3_executable("${f}"
    "${t}"
    USE_PROJECTS SQuIDS
    USE_TOOLS gsl
    WITHOUT_I3_HEADERS)
  set_source_files_properties(${t} COMPILE_FLAGS "-std=c++11")
  add_dependencies(test-bins ${PROJECT_NAME}-${f})
  add_test(NAME "${PROJECT_NAME}::${f}"
           WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
	   COMMAND "${PROJECT_NAME}-${f}")
endforeach()