i3_project(SQuIDS)

i3_add_library(SQuIDS
  private/SQuIDS/*.cpp
  USE_TOOLS gsl
  WITHOUT_I3_HEADERS)
set_source_files_properties(${PROJECT_SOURCE_DIR}/private/${PROJECT_NAME}/const.cpp  COMPILE_FLAGS "-std=c++11")
set_source_files_properties(${PROJECT_SOURCE_DIR}/private/${PROJECT_NAME}/SQUIDS.cpp COMPILE_FLAGS "-std=c++11")
set_source_files_properties(${PROJECT_SOURCE_DIR}/private/${PROJECT_NAME}/SUNalg.cpp COMPILE_FLAGS "-std=c++11 -Warray-bounds")
