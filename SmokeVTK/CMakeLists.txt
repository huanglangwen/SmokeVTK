cmake_minimum_required(VERSION 2.8)
 
PROJECT(WriteVTI)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(WriteVTI MACOSX_BUNDLE WriteVTI)
 
if(VTK_LIBRARIES)
  target_link_libraries(WriteVTI ${VTK_LIBRARIES})
else()
  target_link_libraries(WriteVTI vtkHybrid vtkWidgets)
endif()