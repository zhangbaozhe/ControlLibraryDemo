/*
 * File: acado_invpend_mpc.cpp
 * Description: Use acado MPC on an inverted pendulum.
 *        The position of the inverted pendulum will stay around 0.
 *        The pole of the inverted pendulum will stay still during the control process. 
 *        The control input is force F applied on the cart of the inverted pendulum.
 *        The state can be [x, x_dot, theta, theta_dot]
 * ================================================================================================
 * Author: Baozhe Zhang
 * Date created: Jan, 2023
 */

#include <iostream>
#include <vector>
#include <cmath>

#include <ros/ros.h>
#include <sensor_msgs/JointState.h>
#include <std_msgs/Float64.h>
#include <Eigen/Core>

#include <control_demo/acado_invpend_wrapper.hpp>

using namespace acado_invpend;

const char* CART_NAME = "slider_to_cart";
const char* POLE_NAME = "cart_to_pole";
const char* JOINT_TOPIC = "/invpend/joint_states";
const char* CONTROL_TOPIC = "/invpend/joint1_force_controller/command";

int main(int argc, char **argv)
{
  // Initialize ros node
  ros::init(argc, argv, "mpc_demo_node");
  ros::NodeHandle nh;
  ros::NodeHandle pnh("~");

  // Current estimation
  Eigen::Vector4d CURRENT_EST(Eigen::Vector4d::Zero());
  // Reference
  std::vector<Eigen::Vector4d> REF(kSamples+1);
  for (auto &i : REF) {
    i = Eigen::Vector4d::Zero();
  }

  // Force publisher
  ros::Publisher effort_pub = nh.advertise<std_msgs::Float64>(
    CONTROL_TOPIC, 1
  );

  // Estimation subscriber
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
        CURRENT_EST(1) = msg_ptr->velocity.at(cart_index);
        CURRENT_EST(2) = msg_ptr->position.at(pole_index);
        CURRENT_EST(3) = msg_ptr->velocity.at(pole_index);
      });

  // Wait until connected
  while (!joint_states_sub.getNumPublishers()) {
    ROS_INFO("Waiting ...");
    ros::spinOnce();
  }

  // Define acado mpc controller
  MpcController<double> controller(nh, pnh);

  // Define acado mpc parameter
  MpcParams<double> params;

  ros::Rate freq(200.0);


  while (ros::ok()) {
    // Run the optimization in this iteration
    auto command = controller.run(CURRENT_EST, REF, params);
    // ROS_INFO("EFFORT: %.2f", command.data);

    effort_pub.publish(command);

    freq.sleep();
    ros::spinOnce();
  }


}