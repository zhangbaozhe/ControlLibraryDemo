/*
 * File： ct_quad_circle_mpc.cpp
 * Description: Use `ct_optcon' package to run a demo MPC on iris quadrotor.
 *      The quadrotor will follow a predefined circle trajectory.
 *      The control input is [T, w_x, w_y, w_z].
 *      The state is [x, y, z, x_dot, y_dot, z_dot, q_w, q_x, y_y, q_z].
 * =========================================================================
 * Author: Baozhe Zhang
 * Date created: Dec 12, 2022
 */


/* 
 * This define is necessary since I dont know why control toolbox pacakge cannot find its dependency 
 * `HPIPM' in CMake and then define `HPIPM' during build time. If this is not defined, the bounding
 * box constraints cannot be compiled.
 */

// This defination works for using auto diffrentiate code generating tool provided by control toolbox
#define ADCG_LINEARIZER

#include <chrono>
#include <memory>
#include <iostream>
#include <vector>
#include <cmath>

#include <ct/core/core.h>
#include <ct/core/systems/continuous_time/linear/ADCodegenLinearizer.h>
#include <ct/optcon/optcon.h>
#include <ct/core/integration/Integrator.h>

#include <Eigen/Core>

#include <ros/ros.h>
#include <mavros_msgs/AttitudeTarget.h>
#include <gazebo_msgs/ModelStates.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/Pose.h>

using namespace ct::core;
using namespace ct::optcon;

const size_t state_dim = 10;
const size_t control_dim = 4;
const double g = 9.81;

const size_t SAMPLE_NUM           = 1000;
const size_t STATE_DIM            = 10;
const size_t CONTROL_DIM          = 4;
const double TRAJ_R               = 5;
const double TRAJ_W               = M_PI / 8;
const double TRAJ_DURATION        = 16; // seconds
const double G                    = 9.81;
const double MIN_THRUST           = 2.0;
const double MAX_THRUST           = 20.0;
const double THRUST_PARAM         = 0.7; // thrust when the quadrotor can hover

const char *SETPOINT_ATTITUDE     = "mavros/setpoint_raw/attitude";
const char *GAZEBO_MODELS_TOPIC   = "/gazebo/model_states";
const char *MODEL_NAME            = "iris";



/* ============== CT classes ============== */

// Quadrotor dynamics 
template <typename SCALAR>
class QuadSystem final : public ControlledSystem<state_dim, control_dim, SCALAR>
{
 public: 
  static const size_t STATE_DIM = state_dim;
  static const size_t CONTROL_DIM = control_dim;

  using Base = ControlledSystem<state_dim, control_dim, SCALAR>;
  using time_t = typename Base::time_t;


  // constructor
  QuadSystem(std::shared_ptr<Controller<state_dim, control_dim, SCALAR>> controller = nullptr)
      : Base(controller)
  {
  }

  // copy constructor
  QuadSystem(const QuadSystem &rhs)
      : Base(rhs)
  {
  }

  // deep copy
  QuadSystem *clone() const override { return new QuadSystem(*this); }

  // dynamics
  virtual void computeControlledDynamics(const StateVector<state_dim, SCALAR> &state,
                                 const SCALAR &t,
                                 const ControlVector<control_dim, SCALAR> &control,
                                 StateVector<state_dim, SCALAR> &derivative) override
  {
#ifdef ADCG_LINEARIZER
    derivative(0) = state(3); // vx
    derivative(1) = state(4); // vy
    derivative(2) = state(5); // vz
    derivative(3) = 2 * ( state(6) * state(8) + state(7) * state(9) ) * control(0); // ax
    derivative(4) = 2 * ( state(8) * state(9) - state(6) * state(7) ) * control(0); // ay
    derivative(5) = ( 1 - 2 * state(7) * state(7) - 2 * state(8) * state(8) ) * control(0) - g; // az
    derivative(6) = 0.5 * ( - control(1) * state(7) - control(2) * state(8) - control(3) * state(9));
    derivative(7) = 0.5 * ( control(1) * state(6) + control(3) * state(8) - control(2) * state(9));
    derivative(8) = 0.5 * ( control(2) * state(6) - control(3) * state(7) + control(1) * state(9));
    derivative(9) = 0.5 * ( control(3) * state(6) + control(2) * state(7) - control(1) * state(8));
#else
    Eigen::Quaterniond q_WB(state(6), state(7), state(8), state(9));
    Eigen::Vector3d temp = q_WB * Eigen::Vector3d(0, 0, control(0));
    Eigen::Matrix4d Lambda;
    Lambda << 0, -control(1), -control(2), -control(3),
        control(1), 0, control(3), -control(2),
        control(2), -control(3), 0, control(1),
        control(3), control(2), -control(1), 0;
    derivative(0) = state(3); // vx
    derivative(1) = state(4); // vy
    derivative(2) = state(5); // vz
    derivative(3) = temp(0); // ax
    derivative(4) = temp(1); // ay
    derivative(5) = -9.81 + temp(2); // az
    Eigen::Vector4d temp_vector = 0.5 * Lambda * Eigen::Vector4d(
        q_WB.w(), q_WB.x(), q_WB.y(), q_WB.z());
    derivative(6) = temp_vector(0);
    derivative(7) = temp_vector(1);
    derivative(8) = temp_vector(2);
    derivative(9) = temp_vector(3);
#endif
  }
};

// Quadrotor MPC wrapper
class QuadMPC
{
 public:
  static const size_t STATE_DIM = state_dim;
  static const size_t CONTROL_DIM = control_dim;
  const std::string SOLVER_DIR = "/solver.info";
  const std::string COST_DIR = "/cost.info";

  using OptConProblem_t = ContinuousOptConProblem<state_dim, control_dim>;
  using NLOptConSolver_t = NLOptConSolver<state_dim, control_dim>;
  using State_t = StateVector<state_dim>;
  using Control_t = ControlVector<control_dim>;
  using System_t = QuadSystem<double>; // if use code generation, the scalar needs to be ct::core::ADCodegenLinearizer<state_dim, control_dim>::ADCGScalar
#ifdef ADCG_LINEARIZER
  using SystemADCG_t = QuadSystem<ct::core::ADCGScalar>;
  using CGLinearizer_t = ct::core::ADCodegenLinearizer<state_dim, control_dim>;
#endif
  using Q_t = Eigen::Matrix<double, state_dim, state_dim>;
  using R_t = Eigen::Matrix<double, control_dim, control_dim>;
  using TermQuadratic_t = ct::optcon::TermQuadratic<state_dim, control_dim>;
  using Cost_t = CostFunctionAnalytical<state_dim, control_dim>;
  using Setting_t = NLOptConSettings;
  using BoxConstraint_t = ControlInputConstraint<state_dim, control_dim>;
  using MPC_t = MPC<NLOptConSolver_t>;
  using Integrator_t = Integrator<state_dim>;
  

  Q_t Q_;
  Q_t Q_final_;
  R_t R_;
  State_t x0_;
  Control_t u_max_;
  Control_t u_min_;
  Setting_t settings_;
  Setting_t nloc_settings_mpc_;
  mpc_settings mpc_settings_;
  Time time_horizon_;
  size_t K_;
  double dt_; // time interval
  StateVectorArray<state_dim> state_ref_traj_;
  ControlVectorArray<control_dim> ff_input_;
  std::vector<size_t> inter_term_id_vec_;
  size_t final_term_id_;

  std::shared_ptr<System_t> system_ptr_;
#ifdef ADCG_LINEARIZER
  std::shared_ptr<SystemADCG_t> systemADCG_ptr_;
  std::shared_ptr<CGLinearizer_t> linearizer_ptr_;
#endif
  std::shared_ptr<Cost_t> cost_ptr_;
  std::shared_ptr<OptConProblem_t> opt_prob_ptr_;
  std::shared_ptr<NLOptConSolver_t> solver_ptr_;
  std::shared_ptr<MPC_t> mpc_ptr_;
  std::shared_ptr<Integrator_t> integrator_ptr_;

  // default constructor
  QuadMPC() = delete;

  // constructor
  QuadMPC(const std::string &load_path,
          const Control_t &u_max, const Control_t &u_min, bool verbose = true)
      : u_max_(u_max), u_min_(u_min)
  {
    system_ptr_ = std::make_shared<System_t>();    
#ifdef ADCG_LINEARIZER
    systemADCG_ptr_ = std::make_shared<SystemADCG_t>();
    linearizer_ptr_ = std::make_shared<CGLinearizer_t>(systemADCG_ptr_);
    std::cout << "compiling ..." << std::endl;
    linearizer_ptr_->compileJIT("Quad_ADCGLib");
    std::cout << "... done" << std::endl;
#endif
    integrator_ptr_ = std::make_shared<Integrator_t>(system_ptr_, ct::core::IntegrationType::RK4CT);

    // load optimal control settings
    settings_.load(load_path + SOLVER_DIR, verbose, "nloc");

    loadScalar(load_path + SOLVER_DIR, "time_horizon", time_horizon_);
    loadMatrix(load_path + COST_DIR, "x0", x0_);
    K_ = settings_.computeK(time_horizon_);
    dt_ = time_horizon_ / K_;


    // construct cost function terms
    // construct cost function
    cost_ptr_ = std::make_shared<Cost_t>();
    for (size_t i = 0; i < K_; i++) {
      std::shared_ptr<TermQuadratic_t> term_intermediate(new TermQuadratic_t);
      term_intermediate->loadConfigFile(load_path + COST_DIR, "term0", verbose);
      inter_term_id_vec_.push_back(cost_ptr_->addIntermediateTerm(term_intermediate));
    }
    std::shared_ptr<TermQuadratic_t> term_final(new TermQuadratic_t);
    term_final->loadConfigFile(load_path + COST_DIR, "term1", verbose);
    final_term_id_ = cost_ptr_->addFinalTerm(term_final);

    // construct opt problem
#ifdef ADCG_LINEARIZER
    opt_prob_ptr_ = std::make_shared<OptConProblem_t>(time_horizon_, x0_, system_ptr_, cost_ptr_, linearizer_ptr_);
#else
    opt_prob_ptr_ = std::make_shared<OptConProblem_t>(time_horizon_, x0_, system_ptr_, cost_ptr_);
#endif
    solver_ptr_ = std::make_shared<NLOptConSolver_t>(*opt_prob_ptr_, settings_);

    std::shared_ptr<BoxConstraint_t> control_input_bound(
        new BoxConstraint_t(u_min_, u_max_));
    control_input_bound->setName("ControlInputBound");
    std::shared_ptr<ct::optcon::ConstraintContainerAnalytical<state_dim, control_dim>> input_box_containers(
        new ConstraintContainerAnalytical<state_dim, control_dim>());
    input_box_containers->addIntermediateConstraint(control_input_bound, verbose);
    input_box_containers->initialize();

    State_t init_x = State_t::Zero();
    // hover
    init_x(2) = 1.0;

    StateVectorArray<state_dim> state_ref_traj(K_ + 1, init_x);
    FeedbackArray<state_dim, control_dim> fb_traj(K_, FeedbackMatrix<state_dim, control_dim>::Zero());
    ControlVectorArray<control_dim> ff_traj(K_, ct::core::ControlVector<control_dim>::Zero());

    NLOptConSolver_t::Policy_t init_controller(state_ref_traj, ff_traj, fb_traj, settings_.dt);
    solver_ptr_->setInitialGuess(init_controller);
    solver_ptr_->solve();
    std::cout << "In constructor of the MPC class" << std::endl;;
    std::cout << "Initial solving OK" << std::endl;

    StateFeedbackController<state_dim, control_dim> init_solution = solver_ptr_->getSolution();

    // set constraints
    opt_prob_ptr_->setInputBoxConstraints(input_box_containers);

    // mpc settings
    nloc_settings_mpc_ = solver_ptr_->getSettings();
    nloc_settings_mpc_.max_iterations = 1;
    nloc_settings_mpc_.printSummary = verbose;

    mpc_settings_.stateForwardIntegration_ = false;
    mpc_settings_.postTruncation_ = false;
    mpc_settings_.measureDelay_ = false;
    mpc_settings_.mpc_mode = MPC_MODE::CONSTANT_RECEDING_HORIZON;
    mpc_settings_.coldStart_ = false;

    mpc_ptr_ = std::make_shared<MPC_t>(*opt_prob_ptr_, nloc_settings_mpc_, mpc_settings_);
    mpc_ptr_->setInitialGuess(init_solution);
    system_ptr_->setController(std::shared_ptr<StateFeedbackController<state_dim, control_dim>>(
        new StateFeedbackController<state_dim, control_dim>(init_solution)));
  }

  // run the MPC controller for one iteration
  // and get the predicted trajectory
  bool run(const State_t &estimated_state,
           const StateVectorArray<state_dim> &ref_traj,
           std::shared_ptr<StateFeedbackController<state_dim, control_dim>> new_policy_ptr,
           StateVectorArray<state_dim> &predicted_traj)
  {
    auto start_time = std::chrono::high_resolution_clock::now();

    if (ref_traj.size() != K_ + 1) {
      std::cerr << "Wrong size of reference trajectory" << std::endl;
    }

    auto instance = mpc_ptr_->getSolver().getCostFunctionInstances();
    // update reference
    // there may be multiple instances (multiple threads) see the Issue at GitHub
    for (size_t j = 0; j < instance.size(); j++) {
      for (size_t idx = 0; idx < inter_term_id_vec_.size(); idx++) {
        instance[j]->getIntermediateTermById(
            inter_term_id_vec_.at(idx))->updateReferenceState(ref_traj.at(idx));
      }
      instance[j]->getFinalTermById(final_term_id_)->updateReferenceState(ref_traj.back());
    }

    // ready to solve
    auto current_time = std::chrono::high_resolution_clock::now();
    ct::core::Time t =
        1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(current_time - start_time).count();
    mpc_ptr_->prepareIteration(t);

    ct::core::Time ts_new_policy;
    current_time = std::chrono::high_resolution_clock::now();
    t = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(current_time - start_time).count();

    // copy
    State_t estimated_state_copy = estimated_state;
    // check if the estimated and reference quaternion live in the same hemisphere
    if (estimated_state_copy.segment(6, 4).dot(ref_traj[0].segment(6, 4)) < 0.0) {
      estimated_state_copy.segment(6, 4) = -estimated_state_copy.segment(6, 4);
      ROS_WARN("Quaternion changing");
    }

    bool success = mpc_ptr_->finishIteration(
        estimated_state_copy, t, *new_policy_ptr, ts_new_policy);

    

    // compute the predicted trajectory
    State_t initial_state = estimated_state;
    ct::core::tpl::TimeArray<double> time_array;
    system_ptr_->setController(new_policy_ptr);
    integrator_ptr_->integrate_n_steps(initial_state, t, K_, dt_, predicted_traj, time_array);

    return success;
  }
};


/* ============== CT classes ============== */


/* ============== Utility functions ============== */

StateVectorArray<state_dim>
GenerateCircleTrajectory(size_t sample_num, double R, double omega, double duration, 
    double start_z)
{
  StateVectorArray<state_dim> result;
  Eigen::Matrix3d temp_rotation_WB;
  Eigen::Quaterniond temp_q_WB;
  std::vector<Eigen::Quaterniond> q_vector;
  double t = sqrt(pow(G, 2) + pow(omega * omega * R, 2));
  for (size_t i = 0; i < sample_num; i++) {
    temp_rotation_WB << G / t * cos(omega * (i * duration / sample_num)), 
                        -sin(omega * (i * duration / sample_num)), 
                        -omega * omega * R / t * cos(omega * (i * duration / sample_num)), 
                        G / t * sin(omega * (i * duration / sample_num)), 
                        cos(omega * (i * duration / sample_num)), 
                        -omega * omega * R /t * sin(omega * (i * duration / sample_num)), 
                        omega * omega * R / t, 
                        0, 
                        G / t;
    temp_q_WB = Eigen::Quaterniond(temp_rotation_WB);
    temp_q_WB.normalize();
    q_vector.push_back(temp_q_WB);
    // Quaternion may flip
    if (i != 0) {
      Eigen::Vector4d q_prev(q_vector.at(i-1).w(), q_vector.at(i-1).x(), q_vector.at(i-1).y(), q_vector.at(i-1).z());
      Eigen::Vector4d q_current(q_vector.at(i).w(), q_vector.at(i).x(), q_vector.at(i).y(), q_vector.at(i).z());
      if ((q_prev - q_current).squaredNorm() > (q_prev + q_current).squaredNorm()) {
        q_vector.at(i).w() = -q_vector.at(i).w();
        q_vector.at(i).x() = -q_vector.at(i).x();
        q_vector.at(i).y() = -q_vector.at(i).y();
        q_vector.at(i).z() = -q_vector.at(i).z();
      }
    }
    StateVector<state_dim> point;
    point(0) = R * cos(omega * (i * duration / sample_num));
    point(1) = R * sin(omega * (i * duration / sample_num));
    point(2) = start_z;
    point(3) = -omega * R * sin(omega * (i * duration / sample_num));
    point(4) = omega * R * cos(omega * (i * duration / sample_num));
    point(5) = 0;
    point(6) = q_vector.at(i).w();
    point(7) = q_vector.at(i).x();
    point(8) = q_vector.at(i).y();
    point(9) = q_vector.at(i).z();

    result.push_back(point);
  }
  return result;
}

double compute_distance(
    const StateVector<state_dim> &state_estimate, 
    const StateVector<state_dim> & point)
{
  double r = 0.0;
  r = pow(state_estimate(0) - point(0), 2) + 
      pow(state_estimate(1) - point(1), 2) + 
      pow(state_estimate(2) - point(2), 2);
  return r;
}

StateVectorArray<state_dim>
ExtractReference(size_t window_num, 
    const StateVector<state_dim> &state_estimate, 
    const StateVectorArray<state_dim> &trajectory)
{
  StateVectorArray<state_dim> result;
  size_t target_index = 0;
  auto size = trajectory.size();

  using pair = std::pair<size_t, double>;
  std::vector<pair> distance_pairs;
  distance_pairs.resize(size);

  for (auto i = 0; i < size; i++) {
    distance_pairs.at(i) = 
        pair(i, compute_distance(state_estimate, trajectory.at(i)));
  }

  auto min_iter = std::min_element(distance_pairs.begin(), distance_pairs.end(),  
      [=](const pair &a, const pair &b)
      {
        return a.second < b.second;
      }
  );

  target_index = min_iter->first;

  if (target_index > size - window_num) target_index = 0;

  for (auto i = 0; i < window_num; i++)
    result.push_back(trajectory.at(i + target_index));
  
  return result;
}

nav_msgs::Path
StateVectorArray2Path(const StateVectorArray<state_dim> &trajectory)
{
  nav_msgs::Path result;
  result.header.frame_id = "map";
  result.header.stamp = ros::Time::now();
  for (size_t i = 0; i < trajectory.size(); i++) {
    geometry_msgs::PoseStamped temp_pose;
    temp_pose.header.frame_id = "map";
    temp_pose.header.stamp = ros::Time::now();
    temp_pose.pose.orientation.w = trajectory.at(i)(6);
    temp_pose.pose.orientation.x = trajectory.at(i)(7);
    temp_pose.pose.orientation.y = trajectory.at(i)(8);
    temp_pose.pose.orientation.z = trajectory.at(i)(9);
    temp_pose.pose.position.x = trajectory.at(i)(0);
    temp_pose.pose.position.y = trajectory.at(i)(1);
    temp_pose.pose.position.z = trajectory.at(i)(2);
    result.poses.push_back(temp_pose);
  }
  return result;
}

geometry_msgs::Pose
StateVector2Pose(const StateVector<state_dim> &state)
{
  geometry_msgs::Pose result;
  result.position.x = state(0);
  result.position.y = state(1);
  result.position.z = state(2);
  result.orientation.w = state(6);
  result.orientation.x = state(7);
  result.orientation.y = state(8);
  result.orientation.z = state(9);
  return result;
}

mavros_msgs::AttitudeTarget 
ControlVector2AttitudeTarget(const ControlVector<control_dim> &u)
{
  mavros_msgs::AttitudeTarget msg;
  msg.header.stamp = ros::Time::now();
  msg.type_mask = msg.IGNORE_ATTITUDE;
  msg.thrust = std::min(u(0) * THRUST_PARAM / G, 1.0); // map it to 0 ~ 1
  msg.body_rate.x = u(1);
  msg.body_rate.y = u(2);
  msg.body_rate.z = u(3);
  return msg;
}




/* ============== Utility functions ============== */

int main(int argc, char **argv)
{
  ros::init(argc, argv, "mpc_demo_node");
  ros::NodeHandle nh;

  StateVector<state_dim> CURRENT_QUAD_EST;

  std::string load_path = "/home/vm/control_ws/src/control_demo/parameters";

  ros::Subscriber gazebo_model_states_sub = nh.subscribe<gazebo_msgs::ModelStates>(
      GAZEBO_MODELS_TOPIC, 1,
      [&](const gazebo_msgs::ModelStates::ConstPtr &gz_msg_ptr)
      {
        size_t index = -1;
        for (size_t i = 0; i < gz_msg_ptr->name.size(); i++) {
          if (gz_msg_ptr->name.at(i) == MODEL_NAME)
            index = i; 
        }
        if (index == -1) {
          ROS_ERROR("Cannot find the target");
          ros::shutdown();
          exit(1);
        }
        CURRENT_QUAD_EST(0) = gz_msg_ptr->pose.at(index).position.x;
        CURRENT_QUAD_EST(1) = gz_msg_ptr->pose.at(index).position.y;
        CURRENT_QUAD_EST(2) = gz_msg_ptr->pose.at(index).position.z;
        CURRENT_QUAD_EST(3) = gz_msg_ptr->twist.at(index).linear.x;
        CURRENT_QUAD_EST(4) = gz_msg_ptr->twist.at(index).linear.y;
        CURRENT_QUAD_EST(5) = gz_msg_ptr->twist.at(index).linear.z;
        Eigen::Quaterniond temp_q(gz_msg_ptr->pose.at(index).orientation.w, 
            gz_msg_ptr->pose.at(index).orientation.x, 
            gz_msg_ptr->pose.at(index).orientation.y, 
            gz_msg_ptr->pose.at(index).orientation.z);
        temp_q.normalize();
        CURRENT_QUAD_EST(6) = temp_q.w();
        CURRENT_QUAD_EST(7) = temp_q.x();
        CURRENT_QUAD_EST(8) = temp_q.y();
        CURRENT_QUAD_EST(9) = temp_q.z();
      });

  while (!gazebo_model_states_sub.getNumPublishers()) {
    // 
    ros::spinOnce();
  }
  ros::Publisher mavros_attitude_pub = nh.advertise<mavros_msgs::AttitudeTarget>(
      SETPOINT_ATTITUDE, 1);
  ros::Publisher global_path_pub = nh.advertise<nav_msgs::Path>(
      "demo_mpc/global_path", 1
  );
  ros::Publisher reference_path_pub = nh.advertise<nav_msgs::Path>(
      "demo_mpc/reference_path", 1
  );
  ros::Publisher predicted_path_pub = nh.advertise<nav_msgs::Path>(
      "demo_mpc/predicted_path", 1
  );
  ros::Publisher reference_pose_pub = nh.advertise<geometry_msgs::Pose>(
      "demo_mpc/reference_pose", 1
  );
  // Generate trajectory
  auto trajectory = GenerateCircleTrajectory(
      SAMPLE_NUM, TRAJ_R, TRAJ_W, TRAJ_DURATION, 5
  );

  ControlVector<control_dim> u_max, u_min;
  u_max << MAX_THRUST, M_PI / 2, M_PI / 2, M_PI; 
  u_min << MIN_THRUST, -M_PI / 2, -M_PI / 2, -M_PI; 
  QuadMPC mpc_controller(load_path, u_max, u_min, false);
  size_t window_num = mpc_controller.K_ + 1;

  ros::Rate freq(200.0);


  auto new_policy_ptr = std::make_shared<StateFeedbackController<state_dim, control_dim>>();

  while (ros::ok()) {
    global_path_pub.publish(StateVectorArray2Path(trajectory));

    auto x_ref = ExtractReference(window_num, CURRENT_QUAD_EST, trajectory);
    reference_pose_pub.publish(StateVector2Pose(x_ref[0]));

    reference_path_pub.publish(StateVectorArray2Path(x_ref));

    StateVectorArray<state_dim> x_predicted;
    mpc_controller.run(CURRENT_QUAD_EST, x_ref, new_policy_ptr, x_predicted);
    auto command = new_policy_ptr->uff()[0];

#ifdef ADCG_LINEARIZER
#ifdef ADCG_LOG
    std::cout << "State matrix\n" << mpc_controller.linearizer_ptr_->getDerivativeState(CURRENT_QUAD_EST, command) << std::endl;
    std::cout << "Control matrix\n" << mpc_controller.linearizer_ptr_->getDerivativeControl(CURRENT_QUAD_EST, command) << std::endl;
#endif
#endif

    mavros_attitude_pub.publish(ControlVector2AttitudeTarget(command));
    predicted_path_pub.publish(StateVectorArray2Path(new_policy_ptr->x_ref()));

    freq.sleep();
    ros::spinOnce();
  }

}
