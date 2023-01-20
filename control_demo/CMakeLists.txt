cmake_minimum_required(VERSION 3.0.2)
project(control_demo)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_STANDARD 11)



find_package(catkin REQUIRED COMPONENTS
  ct_core
  ct_optcon
  roscpp
  mavros_msgs
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

set(acado_invpend_mpc_dep_src
  ${PROJECT_SOURCE_DIR}/acado_invpend_model/invpend_mpc_codegen/acado_qpoases_interface.cpp
  ${PROJECT_SOURCE_DIR}/acado_invpend_model/invpend_mpc_codegen/acado_integrator.c
  ${PROJECT_SOURCE_DIR}/acado_invpend_model/invpend_mpc_codegen/acado_solver.c
  ${PROJECT_SOURCE_DIR}/acado_invpend_model/invpend_mpc_codegen/acado_auxiliary_functions.c
  ${PROJECT_SOURCE_DIR}/../acado/external_packages/qpoases/SRC/Bounds.cpp
  ${PROJECT_SOURCE_DIR}/../acado/external_packages/qpoases/SRC/Constraints.cpp
  ${PROJECT_SOURCE_DIR}/../acado/external_packages/qpoases/SRC/CyclingManager.cpp
  ${PROJECT_SOURCE_DIR}/../acado/external_packages/qpoases/SRC/Indexlist.cpp
  ${PROJECT_SOURCE_DIR}/../acado/external_packages/qpoases/SRC/MessageHandling.cpp
  ${PROJECT_SOURCE_DIR}/../acado/external_packages/qpoases/SRC/QProblem.cpp
  ${PROJECT_SOURCE_DIR}/../acado/external_packages/qpoases/SRC/QProblemB.cpp
  ${PROJECT_SOURCE_DIR}/../acado/external_packages/qpoases/SRC/SubjectTo.cpp
  ${PROJECT_SOURCE_DIR}/../acado/external_packages/qpoases/SRC/Utils.cpp
  ${PROJECT_SOURCE_DIR}/../acado/external_packages/qpoases/SRC/EXTRAS/SolutionAnalysis.cpp)

add_library(acado_invpend_mpc_dep STATIC
  ${acado_invpend_mpc_dep_src})
target_include_directories(acado_invpend_mpc_dep PUBLIC
  ${PROJECT_SOURCE_DIR}/acado_invpend_model/invpend_mpc_codegen
  ${PROJECT_SOURCE_DIR}/../acado/external_packages/qpoases
  ${PROJECT_SOURCE_DIR}/../acado/external_packages/qpoases/INCLUDE
  ${PROJECT_SOURCE_DIR}/../acado/external_packages/qpoases/SRC)
target_link_libraries(acado_invpend_mpc_dep ${catkin_LIBRARIES})

# if you have installed mavros, you can compile
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

add_executable(acado_invpend_mpc_node
  src/acado_invpend_mpc.cpp)
target_include_directories(acado_invpend_mpc_node PUBLIC
  ${PROJECT_SOURCE_DIR}/acado_invpend_model)
target_link_libraries(acado_invpend_mpc_node
  acado_invpend_mpc_dep
  ${catkin_LIBRARIES})

