# Code Repo for `Tech Talk: Introducing Two Control Libraries -- Control Toolbox and ACADO Toolkit'

This is the ROS workspace source repo of the tech talk I gave at Huzhou Institude of ZJU in December, 2022.
It contains the following packages: 
- `control-toolbox`: from https://github.com/ethz-adrl/control-toolbox
  - `ct`: meta package
  - `ct_core`: core library
  - `ct_optcon`: optimal control library
- `kindr`: needed by `control-toolbox`
- `invpend_experiment`: containning the necessary description of a inverted pendulum, I forked from https://github.com/linZHank/invpend_experiment 
- `control_demo`: containning demo control programs using the two libraries written in C++
  - `ct_quad_circle_mpc.cpp`: (non-linear) model predictive control for quadrotor tracking a circle using `control-toolbox`. 

## Installation Guide

## Executables 