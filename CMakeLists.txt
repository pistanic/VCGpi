cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(OpenGL/CGAL  LANGUAGES CXX)

#add_definitions(-std=c++0x)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}-std=c++14 -O")
set(CMAKE_PREFIX_PATH /home/nick/SPL)

# Add Frounding to complie options 

add_compile_options("-frounding-math")
# Find OpenGL & GLUT
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(X11 REQUIRED)
# Find CGAL
find_package(CGAL REQUIRED)

#Find SPL
find_package(SPL REQUIRED)


######   Hack might be needed for some systems   #####
if (NOT GLUT_Xmu_LIBRARY OR NOT GLUT_Xmi_LIBRARY)
	if (NOT GLUT_Xmu_LIBRARY)
		set(GLUT_Xmu_LIBRARY "")
	endif()
	if(NOT GLUT_Xmi_LIBRARY)
		set(GLUT_Xmi_LIBRARY "")	
	endif()
	set(GLUT_LIBRARIES "${GLUT_glut_LIBRARIES}")
endif() 
#######################################################

# Include Directories 
include_directories(${GLUT_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIRS})
include_directories(${CGAL_INCLUDE_DIRS})
include_directories(${SPL_INCLUDE_DIRS})
# Link Libraries 

#link_libraries(${GLUT_LIBRARIES} ${OPENGL_LIBRARIES})
link_libraries(glut X11 GL GLU)
link_libraries(${CGAL_LIBRARY} ${GMP_LIBRARIES})
link_libraries(${SPL_LIBRARYS})

# Generate Executable files
#add_executable(triMeshInfo triMeshInfo.cpp)
#add_executable(limitSurfInfo limitSurfInfo.cpp)
add_executable(PItest piGraphics.cpp)
add_executable(chaosLines chaosLines.cpp)
#add_executable(triMeshView triMeshView.cpp)
