cmake_minimum_required(VERSION 3.0.2)
project(control_demo)

add_compile_options(-std=c++11)


find_package(catkin REQUIRED COMPONENTS
  ct_core
  ct_optcon
  roscpp
)

find_package(gflags REQUIRED)

## include blasfeo and hpipm, assumed to be installed in "/opt"
list(APPEND CMAKE_PREFIX_PATH "/opt")
find_package(blasfeo REQUIRED)
find_package(hpipm REQUIRED)
message(STATUS "Found HPIPM and BLASFEO")
set(HPIPM ON)
list(APPEND HPIPM_LIBS hpipm blasfeo)
list(APPEND ct_optcon_COMPILE_DEFINITIONS HPIPM)




catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES control_demo
#  CATKIN_DEPENDS ct_core ct_optcon roscpp
  DEPENDS ct_core ct_optcon gflags
)


include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${blasfeo_INCLUDE_DIRS}
  ${hpipm_INCLUDE_DIRS}
)


add_executable(ct_quad_circle_mpc_node
  src/ct_quad_circle_mpc.cpp)
target_link_libraries(ct_quad_circle_mpc_node
  ${catkin_LIBRARIES} 
  ct_core 
  ct_optcon 
  gflags)

add_executable(ct_invpend_pid_node
  src/ct_invpend_pid.cpp)
target_link_libraries(ct_invpend_pid_node
  ${catkin_LIBRARIES} 
  ct_core 
  ct_optcon)

add_executable(ct_invpend_lqr_node
  src/ct_invpend_lqr.cpp)
target_link_libraries(ct_invpend_lqr_node
  ${catkin_LIBRARIES} 
  ct_core 
  ct_optcon)


add_executable(ct_invpend_mpc_node
  src/ct_invpend_mpc.cpp)
target_link_libraries(ct_invpend_mpc_node
  ${catkin_LIBRARIES} 
  ct_core 
  ct_optcon
  gflags)
