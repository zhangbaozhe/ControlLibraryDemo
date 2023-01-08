# Code Repo for Tech Talk: Introducing Two Control Libraries -- Control Toolbox and ACADO Toolkit

This is the ROS workspace source repo of the tech talk I gave in Jan, 2023.
It contains the following packages: 
- `control-toolbox` (modified from upstream): from https://github.com/ethz-adrl/control-toolbox
  - `ct`: meta package
  - `ct_core`: core library
  - `ct_optcon`: optimal control library
- `kindr`: needed by `control-toolbox`
- `invpend_experiment`: containing the necessary description of a inverted pendulum, I forked from https://github.com/linZHank/invpend_experiment 
- `acado`: containing ACADO source code. This folder will be ignored during catkin compilation 
- `control_demo`: containing demo control programs using the two libraries written in C++
  - `ct_quad_circle_mpc.cpp` -> executable `ct_quad_circle_mpc_node`: (non-linear) model predictive control for quadrotor tracking a circle using `control-toolbox`
  - `ct_invpend_pid.cpp` -> executable `ct_invpend_pid_node`: using one PID feedback controller to control the inverted pendulum stay stable, but not working
  - `ct_invpend_lqr.cpp` -> executable `ct_invpend_lqr_node`: using linear quadratic regulator feedback controller to control the inverted pendulum staying stable, working, but less robust
  - `ct_invpend_mpc`.cpp -> executable `ct_invpend_mpc_node`: (non-linear) model predictive control for the inverted pendulum staying stable, working
  - `acado_invpend_mpc.cpp` -> executable `acado_invpend_mpc_node`: Linear MPC for the inverted pendulum staying stable, working

## Installation Guide
First, make a workspace folder, say `demo_ws`
```bash
$ mkdir demo_ws
```
Next clone this repo as `src` folder
```bash
$ git clone https://github.com/zhangbaozhe/ControlLibraryDemo.git src
```
Then install the following dependencies:
- catkin tools: 
  ```bash
  $ sudo apt install python-catkin-tools
  ```
- Control Toolbox dependencies: 
  - `HPIPM & blasfeo`: Go to `src/control-toolbox/ct`
    ```bash
    $ cd src/control-toolbox/ct
    $ bash install_hpipm.sh
    ```
- invpend_experiment dependencies: 
  ```bash
  $ sudo apt install ros-noetic-ros-control ros-noetic-ros-controllers ros-noetic-gazebo-ros ros-noetic-gazebo-ros-control
  ```

In the workspace folder: 
```bash
$ catkin build -DCMAKE_BUILD_TYPE=Release -j2
```
If you have more than 8GB RAM or even 16GB, you can maybe drop `-j2` option. 
If you want compile the examples and tests in control toolbox you can type
```bash
$ catkin build -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=true -DBUILD_TESTS=true -j2
```
The compilation may take a while, since C++ templates compilation consumes huge resources. 

(Optional) If you want to use the quadrotor demo, you need to use PX4-Autopilot simulation.

## Executables 

Run the inverted pendulum demo: 
In one terminal (in the workspace folder):
```bash
$ source devel/setup.bash
$ roslaunch invpend_control load_invpend.launch
```
In another terminal
```bash
$ source devel/setup.bash
$ rosrun control_demo ct_invpend_mpc_node
```

