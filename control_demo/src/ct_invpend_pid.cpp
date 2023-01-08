/*
 * File: ct_invpend_pid.cpp
 * Description: Use PID class from `ct_core` to run a demo PID controller on an inverted pendulum.
 *        The position of the inverted pendulum will stay around 0.
 *        The pole of the inverted pendulum will staty still during the control process (theta = pi/2). 
 *        The control input is force F applied on the cart of the inverted pendulum.
 *        The state can be [x, x_dot, theta, theta_dot]
 * Note: Only control theta, a SISO controller, not working
 * ================================================================================================
 * Author: Baozhe Zhang
 * Date created: Dec 17, 2022
 */



#include <chrono>
#include <memory>
#include <iostream>
#include <vector>
#include <cmath>

#include <ct/core/core.h>

#include <Eigen/Core>

#include <ros/ros.h>
#include <sensor_msgs/JointState.h>
#include <std_msgs/Float64.h>

using namespace ct::core;

const char* CART_NAME = "slider_to_cart";
const char* POLE_NAME = "cart_to_pole";
const char* JOINT_TOPIC = "/invpend/joint_states";
const char* CONTROL_TOPIC = "/invpend/joint1_force_controller/command";

const size_t state_dim = 4;
const size_t control_dim = 1;
const double G = 9.81;

template <typename SCALAR>
class InvpendSystem final : public ControlledSystem<state_dim, control_dim, SCALAR>
{
 public: 
  static const size_t STATE_DIM = state_dim;
  static const size_t CONTROL_DIM = control_dim;
  using Base = ControlledSystem<state_dim, control_dim, SCALAR>;

  // function alias
  #define SIN sin
  #define COS cos
  #define TAN tan

  // constructor
  InvpendSystem(const SCALAR &cart_mass, 
      const SCALAR &pend_mass, 
      const SCALAR &friction_param, 
      const SCALAR &pend_len, 
      const SCALAR &pend_mass_inertia, 
      std::shared_ptr<Controller<state_dim, control_dim, SCALAR>> controller = nullptr)
      : M_(cart_mass), 
      m_(pend_mass), 
      b_(friction_param), 
      l_(pend_len), 
      I_(pend_mass_inertia), 
      Base(controller)
  {
  }

  // copy constructor
  InvpendSystem(const InvpendSystem &rhs)
      : M_(rhs.M_), 
      m_(rhs.m_), 
      b_(rhs.m_), 
      l_(rhs.l_), 
      I_(rhs.I_), 
      Base(rhs)
  {
  }

  // deep copy
  InvpendSystem *clone() const override { return new InvpendSystem(*this); }

  // dynamics
  virtual void computeControlledDynamics(const StateVector<state_dim, SCALAR> &state, 
      const SCALAR &t, 
      const ControlVector<control_dim, SCALAR> &control, 
      StateVector<state_dim, SCALAR> &derivative) override
  {
    // x dot
    derivative(0) = state(1);
    // theta dot
    derivative(2) = state(3);
    // x ddot
    derivative(1) = 
        1 / ( M_ + m_ - (m_ * m_ * l_ * l_ * COS(state(2)) * COS(state(2))) / (I_ + m_ * l_ * l_) ) * 
        ( control(0) - b_ * state(1) + m_ * l_ * SIN(state(2)) * state(3) * state(3) - m_ * l_ * COS(state(2)) * m_ * G * l_ * SIN(state(2)) / (I_ + m_ * l_ * l_) );
    // theta ddot
    derivative(3) = 1 / (I_ + m_ * l_ * l_) * 
        (m_ * G * l_ * SIN(state(2)) - m_ * l_ * derivative(1) * COS(state(2)));

    // linearized version
    // x ddot
    // derivative(1) = 1 / ( M_ + m_ - m_ * m_ * l_ * l_ / (I_ + m_ * l_ * l_) ) * 
    //     ( control(0) - b_ * state(1) - m_ * m_ * l_ * l_ * G * state(2) / (I_ + m_ * l_ * l_) ); 
    // theta ddot
    // derivative(3) = ( m_ * G * l_ * state(2) - m_ * l_ * derivative(1) ) / (I_ + m_ * l_ * l_);
  }
 private:
  SCALAR M_;
  SCALAR m_;
  SCALAR b_;
  SCALAR l_;
  SCALAR I_;
};



int main(int argc, char **argv)
{

  ros::init(argc, argv, "pid_demo_node");
  ros::NodeHandle nh;

  StateVector<state_dim> CURRENT_EST;

  ros::Publisher effort_pub = nh.advertise<std_msgs::Float64>(
    CONTROL_TOPIC, 1
  );

  ros::Subscriber joint_states_sub = nh.subscribe<sensor_msgs::JointState>(
      JOINT_TOPIC, 1, 
      [&](const sensor_msgs::JointState::ConstPtr &msg_ptr)
      {
        int pole_index = -1;
        int cart_index = -1;
        for (int i = 0; i < msg_ptr->name.size(); i++) {
          if (msg_ptr->name.at(i) == POLE_NAME) 
            pole_index = i;
          if (msg_ptr->name.at(i) == CART_NAME) 
            cart_index = i;
        }
        if (pole_index == -1 || cart_index == -1) {
          ROS_ERROR("Cannot find any target");
          ros::shutdown();
          exit(1);
        }
        CURRENT_EST(0) = msg_ptr->position.at(cart_index);
        CURRENT_EST(1) = msg_ptr->position.at(pole_index);
        CURRENT_EST(2) = msg_ptr->velocity.at(cart_index);
        CURRENT_EST(3) = msg_ptr->velocity.at(pole_index);
      });

  while (!joint_states_sub.getNumPublishers()) {
    ros::spinOnce();
  }

  PIDController<double>::parameters_t pid_params;
  pid_params.dt = 0.01;
  pid_params.k_p = 10;
  pid_params.k_i = 1;
  pid_params.k_d = 1;
  pid_params.uMax = 2;
  pid_params.uMin = -2;
  PIDController<double>::setpoint_t pid_setpoint;
  pid_setpoint.stateDesired_ = 0;
  PIDController<double> pid_theta_controller(pid_params, pid_setpoint);

  ros::Rate freq(100.0);
  auto start = ros::Time::now().toSec();
  while (ros::ok()) {
    ROS_INFO("EST: angle: %.2f", CURRENT_EST(1));
    auto now = ros::Time::now().toSec() - start;
    auto control = pid_theta_controller.computeControl(CURRENT_EST(1), now);
    std_msgs::Float64 msg;
    msg.data = control;
    effort_pub.publish(msg);

    freq.sleep();
    ros::spinOnce();
  }
  std_msgs::Float64 msg;
  msg.data = 0;
  effort_pub.publish(msg);

  return 0;

}