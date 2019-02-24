cmake_minimum_required(VERSION 2.8)
project(KT C)

set(GKLIB_PATH "GKlib" CACHE PATH "path to GKlib")
set(SHARED FALSE CACHE BOOL "build a shared library")

if(MSVC)
  set(KT_INSTALL FALSE)
else()
  set(KT_INSTALL TRUE)
endif()

# Add GKlib's cmake settings
include(${GKLIB_PATH}/GKlibSystem.cmake)

# Source files
file(GLOB kt_sources *.c)

# Include directories
include_directories(${GKLIB_PATH})
include_directories(.)
include_directories(~/local/include)

# Link directories 
link_directories(~/local/lib)

# Build
add_executable(kt ${GKlib_sources} ${kt_sources})

foreach(prog kt)
  target_link_libraries(${prog} m)
endforeach(prog)

# Install
if(KT_INSTALL)
  install(TARGETS kt 
    RUNTIME DESTINATION bin)
endif()

# Try to find subversion revision.
set(SVNREV "")
file(TO_NATIVE_PATH ${PROJECT_SOURCE_DIR}/.svn svn_dir)
if(IS_DIRECTORY ${svn_dir})
  include(FindSubversion)
  if(Subversion_FOUND)
    Subversion_WC_INFO(${PROJECT_SOURCE_DIR} kt)
    set(SVNREV ${kt_WC_REVISION})
  endif(Subversion_FOUND)
endif()
add_definitions(-DSVNINFO="${SVNREV}")
