/*
 * File: acado_invpend_wrapper.hpp
 * Description: This all-in-one header file provides all the necessary 
 *              interfaces with generated ACADO C code.
 *              This program is inspired by the code in rpg_mpc.
 * ================================================================================================
 * Author: Baozhe Zhang
 * Date created: Dec 25, 2022
 */

#pragma once

#include <vector>
#include <thread>
#include <memory>

#include <Eigen/Eigen>
#include <ros/ros.h>
#include "std_msgs/Float64.h"

namespace acado_invpend
{

#include "invpend_mpc_codegen/acado_auxiliary_functions.h"
#include "invpend_mpc_codegen/acado_common.h"

static constexpr int kSamples = ACADO_N;      // number of samples
static constexpr int kStateSize = ACADO_NX;   // number of states
static constexpr int kRefSize = ACADO_NY;     // number of reference states
static constexpr int kEndRefSize = ACADO_NYN; // number of end reference states
static constexpr int kInputSize = ACADO_NU;   // number of inputs
static constexpr int kCostSize = ACADO_NY - ACADO_NU; // number of state costs

ACADOvariables acadoVariables;
ACADOworkspace acadoWorkspace;

enum class STATE {
  kX = 0, kXDot, kTheta, kThetaDot
};

enum class INPUT {
  kU = 0
};

// Parameter utility functions
template <typename T>
bool getParam(const std::string &name, T &parameter, 
    const ros::NodeHandle &pnh = ros::NodeHandle("~"))
{
  if (pnh.getParam(name, parameter)) {
    ROS_INFO_STREAM(
        "[" << pnh.getNamespace() << "]" << name << " = " << parameter);
    return true;
  }
  ROS_ERROR_STREAM(
      "[" << pnh.getNamespace() << "] Could not load parameter "
      << pnh.getNamespace() << "/"<< name);
  return false;
}

template<typename T>
bool getParam(const std::string& name, T& parameter, const T& defaultValue,
              const ros::NodeHandle& pnh = ros::NodeHandle("~"))
{
  if (pnh.getParam(name, parameter))
    ROS_INFO_STREAM(
        "[" << pnh.getNamespace() << "] " << name << " = " << parameter);
  else {
    parameter = defaultValue;
    ROS_WARN_STREAM(
        "[" << pnh.getNamespace() << "] " << name << " = " << parameter
        << " set to default value");
  }
  return true;
}

// Declarations

template <typename T>
class MpcParams
{
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  MpcParams() : 
    changed_(false), 
    print_info_(false), 
    u_min_(0.0), 
    u_max_(0.0), 
    Q_(Eigen::Matrix<T, kCostSize, kCostSize>::Zero()), 
    R_(Eigen::Matrix<T, kInputSize, kInputSize>::Zero())
  {
  }

  ~MpcParams() {}

  bool loadParameters(ros::NodeHandle &pnh)
  {
    // utility macros
#define GET_PARAM(name) if(!getParam(#name, name, pnh)) return false
#define GET_PARAM_(name) if (!getParam(#name, name##_, pnh)) return false

    T Q_x, Q_x_dot, Q_theta, Q_theta_dot;
    GET_PARAM(Q_x);
    GET_PARAM(Q_x_dot);
    GET_PARAM(Q_theta);
    GET_PARAM(Q_theta_dot);

    if (Q_x <= 0.0 || Q_x_dot < 0.0 || Q_theta <= 0.0 || Q_theta_dot < 0.0) {
      ROS_ERROR("MPC: State cost Q has negative entries");
      return false;
    }
    T R_u;
    GET_PARAM(R_u);
    if (R_u <= 0.0) {
      ROS_ERROR("MPC: Input cost R has negative entries!");
      return false;
    }

    Q_ = (Eigen::Matrix<T, kCostSize, 1>() << 
        Q_x, Q_x_dot, Q_theta, Q_theta_dot).finished().asDiagonal();
    R_ = (Eigen::Matrix<T, kInputSize, 1>() << 
        R_u).finished().asDiagonal();

    GET_PARAM_(u_min);
    GET_PARAM_(u_max);
    if (u_min_ >= u_max_ || u_min_ > 0.0 || u_max_ < 0.0) {
      ROS_ERROR("MPC: Limits not correct");
      return false;
    }

    getParam("print_info", print_info_, false, pnh);
    if (print_info_) ROS_INFO("MPC: Informative printing enabled.");

    changed_ = true;
#undef GET_PARAM
#undef GET_PARAM_
    return true;
  }

  bool changed_;
  bool print_info_;
  T u_min_;
  T u_max_;
  Eigen::Matrix<T, kCostSize, kCostSize> Q_;
  Eigen::Matrix<T, kInputSize, kInputSize> R_;

};


template <typename T>
class MpcWrapper
{
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  MpcWrapper();

  MpcWrapper(
    const Eigen::Ref<const Eigen::Matrix<T, kCostSize, kCostSize>> Q,
    const Eigen::Ref<const Eigen::Matrix<T, kInputSize, kInputSize>> R);

  bool setCosts(
    const Eigen::Ref<const Eigen::Matrix<T, kCostSize, kCostSize>> Q,
    const Eigen::Ref<const Eigen::Matrix<T, kInputSize, kInputSize>> R);

  bool setLimits(T u_min, T u_max);

  bool setReferencePose(
    const Eigen::Ref<const Eigen::Matrix<T, kStateSize, 1>> state);

  bool setTrajectory(
    const Eigen::Ref<const Eigen::Matrix<T, kStateSize, kSamples+1>> states,
    const Eigen::Ref<const Eigen::Matrix<T, kInputSize, kSamples+1>> inputs);

  bool solve(const Eigen::Ref<const Eigen::Matrix<T, kStateSize, 1>> state);

  bool update(const Eigen::Ref<const Eigen::Matrix<T, kStateSize, 1>> state,
              bool do_preparation = true);

  bool prepare();

  void getState(const int node_index,
    Eigen::Ref<Eigen::Matrix<T, kStateSize, 1>> return_state);

  void getStates(
    Eigen::Ref<Eigen::Matrix<T, kStateSize, kSamples+1>> return_states);

  void getInput(const int node_index,
    Eigen::Ref<Eigen::Matrix<T, kInputSize, 1>> return_input);

  void getInputs(
    Eigen::Ref<Eigen::Matrix<T, kInputSize, kSamples>> return_input);

  T getTimestep() { return dt_; }
 private: 
  Eigen::Map<Eigen::Matrix<float, kRefSize, kSamples, Eigen::ColMajor>>
    acado_reference_states_{acadoVariables.y};

  Eigen::Map<Eigen::Matrix<float, kEndRefSize, 1, Eigen::ColMajor>>
    acado_reference_end_state_{acadoVariables.yN};

  Eigen::Map<Eigen::Matrix<float, kStateSize, 1, Eigen::ColMajor>>
    acado_initial_state_{acadoVariables.x0};

  Eigen::Map<Eigen::Matrix<float, kStateSize, kSamples+1, Eigen::ColMajor>>
    acado_states_{acadoVariables.x};

  Eigen::Map<Eigen::Matrix<float, kInputSize, kSamples>>
    acado_inputs_{acadoVariables.u};

  Eigen::Map<Eigen::Matrix<float, kRefSize, kRefSize * kSamples>>
    acado_W_{acadoVariables.W};

  Eigen::Map<Eigen::Matrix<float, kEndRefSize, kEndRefSize>>
    acado_W_end_{acadoVariables.WN};

  Eigen::Map<Eigen::Matrix<float, 4, kSamples, Eigen::ColMajor>>
    acado_lower_bounds_{acadoVariables.lbValues};

  Eigen::Map<Eigen::Matrix<float, 4, kSamples, Eigen::ColMajor>>
    acado_upper_bounds_{acadoVariables.ubValues};
  
  Eigen::Map<Eigen::Matrix<float, 1, kSamples>>
    acado_state_lower_bounds_{acadoVariables.lbValues};

  Eigen::Map<Eigen::Matrix<float, 1, kSamples>>
    acado_state_upper_bounds_{acadoVariables.ubValues};

  Eigen::Matrix<T, kRefSize, kRefSize> W_ = 
      (Eigen::Matrix<T, kRefSize, 1>() << 100, 1, 50000, 1, 1).finished().asDiagonal();

  Eigen::Matrix<T, kEndRefSize, kEndRefSize> WN_ = 
      W_.block(0, 0, kEndRefSize, kEndRefSize);
  
  bool acado_is_prepared_{false};

  const T dt_{0.01};

  const Eigen::Matrix<real_t, kInputSize, 1> kInitialInput_ = 
      (Eigen::Matrix<real_t, kInputSize, 1>() << 0.0).finished();
};

template <typename T>
class MpcController
{
 public: 
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  static_assert(kStateSize == 4, 
      "MpcController: Wrong model size. Number of states does not match.");
  static_assert(kInputSize == 1, 
      "MpcController: Wrong model size. Number of inputs does not match.");

  using Estimate_t = Eigen::Matrix<T, kStateSize, 1>;
  using Trajectory_t = std::vector<Eigen::Matrix<T, kStateSize, 1> >;
  using Control_t = std_msgs::Float64;

  MpcController(const ros::NodeHandle &nh, const ros::NodeHandle &pnh);

  MpcController() : MpcController(ros::NodeHandle(), ros::NodeHandle("~")) {};

  Control_t off();

  Control_t run(const Estimate_t &state_estimate, 
      const Trajectory_t &reference_trajectory, 
      const MpcParams<T> &params);

 private:

  Control_t updateControlCommand(const Estimate_t &state, const Control_t &input, 
      ros::Time &time);

  void preparationThread();

  bool setNewParams(const MpcParams<T> &params);

  ros::NodeHandle nh_;
  ros::NodeHandle pnh_;

  MpcParams<T> params_;

  MpcWrapper<T> mpc_wrapper_;

  std::thread preparation_thread_;

  T timing_feedback_ = 0, timing_preparation_ = 0;
  bool solve_from_scratch_;
  Eigen::Matrix<T, kStateSize, 1> est_state_;
  Eigen::Matrix<T, kStateSize, kSamples + 1> reference_states_;
  Eigen::Matrix<T, kInputSize, kSamples + 1> reference_inputs_;
  Eigen::Matrix<T, kStateSize, kSamples + 1> predicted_states_;
  Eigen::Matrix<T, kInputSize, kSamples> predicted_inputs_;
};





// Implementations

template <typename T>
MpcWrapper<T>::MpcWrapper()
{
  memset(&acadoWorkspace, 0, sizeof(acadoWorkspace));
  memset(&acadoVariables, 0, sizeof(acadoVariables));

  acado_initializeSolver();

  const Eigen::Matrix<T, kStateSize, 1> init_state = 
      (Eigen::Matrix<T, kStateSize, 1>() << 0.0, 0.0, 0.0, 0.0).finished();

  acado_initial_state_ = init_state.template cast<float>();

  acado_states_ = init_state.replicate(1, kSamples + 1).template cast<float>();
  acado_inputs_ = kInitialInput_.replicate(1, kSamples).template cast<float>();

  acado_reference_states_.block(0, 0, kStateSize, kSamples) = 
      init_state.replicate(1, kSamples).template cast<float>();

  acado_reference_states_.block(kStateSize, 0, kCostSize-kStateSize, kSamples) =
      Eigen::Matrix<float, kCostSize-kStateSize, kSamples>::Zero();
  
  // acado_reference_states_.block(kCostSize, 0, kInputSize, kSamples) = 
  //     kInitialInput_.replicate(1, kSamples);
  acado_reference_states_.block(kCostSize, 0, kInputSize, kSamples) = 
      Eigen::Matrix<float, kInputSize, kSamples>::Zero();
  
  acado_reference_end_state_.segment(0, kStateSize) = 
      init_state.template cast<float>();
  
  acado_reference_end_state_.segment(kStateSize, kCostSize-kStateSize) = 
      Eigen::Matrix<float, kCostSize-kStateSize, 1>::Zero();


  if (!acado_W_.trace() > 0.0) {
    acado_W_ = W_.replicate(1, kSamples).template cast<float>();
    acado_W_end_ = WN_.template cast<float>();
  }

  acado_initializeNodesByForwardSimulation();
  acado_preparationStep();
  acado_is_prepared_ = true;
}

template <typename T>
MpcWrapper<T>::MpcWrapper(
    const Eigen::Ref<const Eigen::Matrix<T, kCostSize, kCostSize>> Q,
    const Eigen::Ref<const Eigen::Matrix<T, kInputSize, kInputSize>> R) :
    MpcWrapper()
{
  setCosts(Q, R);
}

template <typename T>
bool MpcWrapper<T>::setCosts(
    const Eigen::Ref<const Eigen::Matrix<T, kCostSize, kCostSize>> Q,
    const Eigen::Ref<const Eigen::Matrix<T, kInputSize, kInputSize>> R)
{
  
  W_.block(0, 0, kCostSize, kCostSize) = Q;
  W_.block(kCostSize, kCostSize, kInputSize, kInputSize) = R;
  WN_ = W_.block(0, 0, kCostSize, kCostSize); 
  for(int i = 0; i < kSamples; i++)
  { 
    acado_W_.block(0, i*kRefSize, kCostSize, kCostSize) =
      W_.block(0, 0, kCostSize, kCostSize).template cast<float>();
    acado_W_.block(kCostSize, i*kRefSize+kCostSize, kInputSize, kInputSize) =
      W_.block(kCostSize, kCostSize, kInputSize, kInputSize).template cast<float>();
  } 
  acado_W_end_ = WN_.template cast<float>();
  std::cout << acado_W_ << std::endl;
  return true;
}

template <typename T>
bool MpcWrapper<T>::setLimits(T u_min, T u_max)
{
  if (u_min >= 0.0 || u_min > u_max) {
    ROS_ERROR("MPC: Minimal effort is not set properly, not changed.");
    return false;
  }

  if (u_max <= 0.0 || u_min > u_max) {
    ROS_ERROR("MPC: Maximal effort is not set properly, not changed.");
    return false;
  }


  acado_lower_bounds_ = Eigen::Matrix<float, kStateSize, kSamples>::Ones() * u_min;

  acado_upper_bounds_ = Eigen::Matrix<float, kStateSize, kSamples>::Ones() * u_max;
  
  return true;
}

template <typename T>
bool MpcWrapper<T>::setReferencePose(
    const Eigen::Ref<const Eigen::Matrix<T, kStateSize, 1>> state)
{
  acado_reference_states_.block(0, 0, kStateSize, kSamples) = 
      state.replicate(1, kSamples).template cast<float>();
  
  acado_reference_states_.block(kStateSize, 0, kCostSize-kStateSize, kSamples) = 
      Eigen::Matrix<float, kCostSize-kStateSize, kSamples>::Zero();

  acado_reference_states_.block(kCostSize, 0, kInputSize, kSamples) = 
      kInitialInput_.replicate(1, kSamples);
  
  acado_reference_end_state_.segment(0, kStateSize) = 
      state.template cast<float>();
  
  acado_reference_end_state_.segment(kStateSize, kCostSize-kStateSize) = 
      Eigen::Matrix<float, kCostSize-kStateSize, 1>::Zero();
  
  acado_initializeNodesByForwardSimulation();
  return true;
}

template <typename T>
bool MpcWrapper<T>::setTrajectory(
    const Eigen::Ref<const Eigen::Matrix<T, kStateSize, kSamples+1>> states, 
    const Eigen::Ref<const Eigen::Matrix<T, kInputSize, kSamples+1>> inputs)
{
  Eigen::Map<Eigen::Matrix<float, kRefSize, kSamples, Eigen::ColMajor>>
      y(const_cast<float*>(acadoVariables.y));
  
  acado_reference_states_.block(0, 0, kStateSize, kSamples) =
      states.block(0, 0, kStateSize, kSamples).template cast<float>();

  acado_reference_states_.block(kStateSize, 0, kCostSize-kStateSize, kSamples) =
      Eigen::Matrix<float, kCostSize-kStateSize, kSamples>::Zero();

  acado_reference_states_.block(kCostSize, 0, kInputSize, kSamples) =
    inputs.block(0, 0, kInputSize, kSamples).template cast<float>();

  acado_reference_end_state_.segment(0, kStateSize) =
      states.col(kSamples).template cast<float>();

  return true;
}

template <typename T>
bool MpcWrapper<T>::solve(
    const Eigen::Ref<const Eigen::Matrix<T, kStateSize, 1>> state)
{
  acado_states_ = state.replicate(1, kSamples+1).template cast<float>();

  acado_inputs_ = kInitialInput_.replicate(1, kSamples);

  return update(state);
}

template <typename T>
bool MpcWrapper<T>::update(
    const Eigen::Ref<const Eigen::Matrix<T, kStateSize, 1>> state,
    bool do_preparation)
{
  if (!acado_is_prepared_) {
    ROS_WARN("MPC: Solver was triggered without preparation, abort!");
    return false;
  }

  acado_initial_state_ = state.template cast<float>();

  // Perform feedback step and reset preparation check.
  acado_feedbackStep();
  acado_is_prepared_ = false;

  // Prepare if the solver if wanted
  if (do_preparation) {
    acado_preparationStep();
    acado_is_prepared_ = true;
  }

  return true;
}

template <typename T>
bool MpcWrapper<T>::prepare()
{
  acado_preparationStep();
  acado_is_prepared_ = true;

  return true;
}

template <typename T>
void MpcWrapper<T>::getState(const int node_index,
    Eigen::Ref<Eigen::Matrix<T, kStateSize, 1>> return_state)
{
  return_state = acado_states_.col(node_index).cast<T>();
}

template <typename T>
void MpcWrapper<T>::getStates(
    Eigen::Ref<Eigen::Matrix<T, kStateSize, kSamples+1>> return_states)
{
  return_states = acado_states_.cast<T>();
}

template <typename T>
void MpcWrapper<T>::getInput(const int node_index,
    Eigen::Ref<Eigen::Matrix<T, kInputSize, 1>> return_input)
{
  return_input = acado_inputs_.col(node_index).cast<T>();
}

template <typename T>
void MpcWrapper<T>::getInputs(
    Eigen::Ref<Eigen::Matrix<T, kInputSize, kSamples>> return_inputs)
{
  return_inputs = acado_inputs_.cast<T>();
}

template <typename T>
MpcController<T>::MpcController(
    const ros::NodeHandle &nh, const ros::NodeHandle &pnh) : 
    nh_(nh), pnh_(pnh)
{
  if (!params_.loadParameters(pnh_)) {
    ROS_ERROR("[%s] Could not load parameters.", pnh_.getNamespace().c_str());
    ros::shutdown();
    return;
  }
  setNewParams(params_);

  solve_from_scratch_ = true;
  preparation_thread_ = std::thread(&MpcWrapper<T>::prepare, mpc_wrapper_);
}

template <typename T>
typename MpcController<T>::Control_t MpcController<T>::off()
{
  Control_t command;
  command.data = 0.0;
  return command;
}

template <typename T>
typename MpcController<T>::Control_t MpcController<T>::run(
    const Estimate_t &state_estimate, 
    const Trajectory_t &reference_trajectory, 
    const MpcParams<T> &params)
{
  ros::Time call_time = ros::Time::now();
  const clock_t start = clock();
  if (params.changed_) {
    params_ = params;
    setNewParams(params_);
  }

  preparation_thread_.join();

  est_state_ = state_estimate;
  for (size_t i = 0; i < kSamples+1; i++) {
    reference_states_.col(i) = reference_trajectory.at(i);
  }
  reference_inputs_.setZero();

  static const bool do_preparation_step(false);

  mpc_wrapper_.setTrajectory(reference_states_, reference_inputs_);
  if (solve_from_scratch_) {
    ROS_INFO("Solving MPC with initial zero input as guess");
    mpc_wrapper_.solve(est_state_);
    solve_from_scratch_ = false;
  } else {
    mpc_wrapper_.update(est_state_, do_preparation_step);
  }
  mpc_wrapper_.getStates(predicted_states_);
  mpc_wrapper_.getInputs(predicted_inputs_);
  std::cout << acado_getObjective() << std::endl;

  preparation_thread_ = std::thread(&MpcController<T>::preparationThread, this);
  const clock_t end = clock();
  timing_feedback_ = 0.9 * timing_feedback_ +
                     0.1 * double(end - start) / CLOCKS_PER_SEC;
  if (params_.print_info_)
    ROS_INFO_THROTTLE(1.0, "MPC Timing: Latency: %1.1f ms  |  Total: %1.1f ms",
                      timing_feedback_ * 1000, (timing_feedback_ + timing_preparation_) * 1000);

  Control_t temp;
  temp.data = predicted_inputs_.col(0)(0);
  // Return the input control command.
  return updateControlCommand(predicted_states_.col(0),
                              temp,
                              call_time);
}

template <typename T>
typename MpcController<T>::Control_t MpcController<T>::updateControlCommand(
    const Estimate_t &state, const Control_t &input, 
    ros::Time &time)
{
  T input_bounded = (T)input.data;

  input_bounded = std::max(params_.u_min_, 
      std::min(params_.u_max_, input_bounded));
  
  Control_t command;
  command.data = input_bounded;
  return command;
}

template <typename T>
void MpcController<T>::preparationThread()
{
  const clock_t start = clock();

  mpc_wrapper_.prepare();

  // Timing
  const clock_t end = clock();
  timing_preparation_ = 0.9 * timing_preparation_ +
                        0.1 * double(end - start) / CLOCKS_PER_SEC;
}

template <typename T>
bool MpcController<T>::setNewParams(const MpcParams<T> &params)
{
  mpc_wrapper_.setCosts(params.Q_, params.R_);
  mpc_wrapper_.setLimits(params.u_min_, params.u_max_);
  return true;
}

} // namespace acado_invpend
