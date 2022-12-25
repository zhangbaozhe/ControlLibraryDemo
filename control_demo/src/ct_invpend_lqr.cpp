/*
 * File: ct_invpend_lqr.cpp
 * Description: on an inverted pendulum.
 *        The position of the inverted pendulum will stay around 0.
 *        The pole of the inverted pendulum will staty still during the control process. 
 *        The control input is force F applied on the cart of the inverted pendulum.
 *        The state can be [x, x_dot, theta, theta_dot]
 * ================================================================================================
 * Author: Baozhe Zhang
 * Date created: Dec 18, 2022
 */



#include <chrono>
#include <memory>
#include <iostream>
#include <vector>
#include <cmath>

#include <ct/core/core.h>
#include <ct/optcon/optcon.h>

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

double constrainAngle(double x){
    x = fmod(x + M_PI, 2*M_PI);
    if (x < 0)
        x += 2*M_PI;
    return x - M_PI;
}


int main(int argc, char **argv)
{
  const double M(0.5);
  const double m(0.2);
  const double b(0.1);
  const double l(0.3);
  const double I(0.006);

  ros::init(argc, argv, "lqr_demo_node");
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
        CURRENT_EST(1) = msg_ptr->velocity.at(cart_index);
        CURRENT_EST(2) = msg_ptr->position.at(pole_index);
        CURRENT_EST(3) = msg_ptr->velocity.at(pole_index);
      });

  while (!joint_states_sub.getNumPublishers()) {
    ros::spinOnce();
  }

  std::shared_ptr<InvpendSystem<double>> invpend_sys_ptr(
      new InvpendSystem<double>(M, m, b, l, I));


  ct::core::SystemLinearizer<state_dim, control_dim> Linearizer(invpend_sys_ptr);

  StateVector<state_dim> x;
  x.setZero();
  ControlVector<control_dim> u;
  u.setZero();

  double t = 0.0;
  auto A = Linearizer.getDerivativeState(x, u, t);
  auto B = Linearizer.getDerivativeControl(x, u, t);


  ct::optcon::TermQuadratic<state_dim, control_dim> quadraticCost;
  Eigen::Matrix<double, 4, 4> Q;
  Q << 100.0, 0, 0, 0, 
        0, 1.0, 0, 0, 
        0, 0, 50000.0, 0, 
        0, 0, 0, 1;
  Eigen::Matrix<double, 1, 1> R;
  R << 1.0;
  ct::optcon::LQR<state_dim, control_dim> lqrSolver;
  ct::core::FeedbackMatrix<state_dim, control_dim> K;

  std::cout << "A: " << std::endl << A << std::endl << std::endl;
  std::cout << "B: " << std::endl << B << std::endl << std::endl;
  std::cout << "Q: " << std::endl << Q << std::endl << std::endl;
  std::cout << "R: " << std::endl << R << std::endl << std::endl;
  lqrSolver.compute(Q, R, A, B, K);
  std::cout << "LQR gain matrix:" << std::endl << K << std::endl;
  // K(0, 2) = -K(0, 2);
  // K(0, 3) = -K(0, 3);
  
  
  ros::Rate freq(100.0);
  auto start = ros::Time::now().toSec();
  while (ros::ok()) {

    u = -K * CURRENT_EST;

    u(0) = std::max(-2.0, u(0));
    u(0) = std::min(2.0, u(0));

    std_msgs::Float64 msg;
    msg.data = u(0);
    effort_pub.publish(msg);

    freq.sleep();
    ros::spinOnce();
  }


  return 0;

}