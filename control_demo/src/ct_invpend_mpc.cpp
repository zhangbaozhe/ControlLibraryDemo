/*
 * File: ct_invpend_mpc.cpp
 * Description: on an inverted pendulum.
 *        The position of the inverted pendulum will stay around 0.
 *        The pole of the inverted pendulum will staty still during the control process. 
 *        The control input is force F applied on the cart of the inverted pendulum.
 *        The state can be [x, x_dot, theta, theta_dot]
 * ================================================================================================
 * Author: Baozhe Zhang
 * Date created: Dec 18, 2022
 */



#define ADCG_LINEARIZER
#undef ADCG_LINEARIZER

#include <chrono>
#include <memory>
#include <iostream>
#include <vector>
#include <cmath>

#include <ct/core/core.h>
#include <ct/optcon/optcon.h>
#include <cppad/cppad.hpp>

#include <Eigen/Core>

#include <ros/ros.h>
#include <sensor_msgs/JointState.h>
#include <std_msgs/Float64.h>

using namespace ct::core;
using namespace ct::optcon;

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
#ifdef ADCG_LINEARIZER
  #define SIN CppAD::sin
  #define COS CppAD::cos
  #define TAN CppAD::tan
#else
  #define SIN sin
  #define COS cos
  #define TAN tan
#endif

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
    std::cout << "I'm copied" << std::endl;
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

class InvpendMPC
{
 public: 
  static const size_t STATE_DIM = state_dim;
  static const size_t CONTROL_DIM = control_dim;
  const std::string SOLVER_DIR = "/invpend_solver.info";
  const std::string COST_DIR = "/invpend_cost.info";

  using OptConProblem_t = ContinuousOptConProblem<state_dim, control_dim>;
  using NLOptConSolver_t = NLOptConSolver<state_dim, control_dim>;
  using State_t = StateVector<state_dim>;
  using Control_t = ControlVector<control_dim>;
  using System_t = InvpendSystem<double>; // if use code generation, the scalar needs to be ct::core::ADCodegenLinearizer<state_dim, control_dim>::ADCGScalar
#ifdef ADCG_LINEARIZER
  using SystemADCG_t = InvpendSystem<ct::core::ADCGScalar>;
  using CGLinearizer_t = ct::core::ADCodegenLinearizer<state_dim, control_dim>;
#endif
  using Q_t = Eigen::Matrix<double, state_dim, state_dim>;
  using R_t = Eigen::Matrix<double, control_dim, control_dim>;
  using TermQuadratic_t = ct::optcon::TermQuadratic<state_dim, control_dim>;
  using Cost_t = CostFunctionAnalytical<state_dim, control_dim>;
  using Setting_t = NLOptConSettings;
  using BoxConstraint_t = ControlInputConstraint<state_dim, control_dim>;
  using MPC_t = MPC<NLOptConSolver_t>;
  

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

  // default constructor
  InvpendMPC() = delete;

  // constructor
  InvpendMPC(const std::string &load_path,
          const Control_t &u_max, const Control_t &u_min, 
          const std::shared_ptr<System_t> &sys_double, 
#ifdef ADCG_LINEARIZER
          const std::shared_ptr<SystemADCG_t> &sys_adcg, 
#endif
          bool verbose = false)
      : u_max_(u_max), u_min_(u_min) 
  {
    system_ptr_ = sys_double;    
#ifdef ADCG_LINEARIZER
    systemADCG_ptr_ = sys_adcg;
    linearizer_ptr_ = std::make_shared<CGLinearizer_t>(systemADCG_ptr_);
    std::cout << "compiling ..." << std::endl;
    linearizer_ptr_->compileJIT("invpend_ADCGLib");
    std::cout << "... done" << std::endl;
#endif

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

    std::shared_ptr<BoxConstraint_t> control_input_bound(
        new BoxConstraint_t(u_min_, u_max_));
    control_input_bound->setName("ControlInputBound");
    std::shared_ptr<ct::optcon::ConstraintContainerAnalytical<state_dim, control_dim>> input_box_containers(
        new ConstraintContainerAnalytical<state_dim, control_dim>());
    input_box_containers->addIntermediateConstraint(control_input_bound, verbose);
    input_box_containers->initialize();
    // set constraints
    opt_prob_ptr_->setInputBoxConstraints(input_box_containers);

    solver_ptr_ = std::make_shared<NLOptConSolver_t>(*opt_prob_ptr_, settings_);

    State_t init_x = State_t::Zero();
    init_x(0) = 0.0;
    // init_x(2) = -3.14 / 2.0;

    StateVectorArray<state_dim> state_ref_traj(K_ + 1, init_x);
    FeedbackArray<state_dim, control_dim> fb_traj(K_, FeedbackMatrix<state_dim, control_dim>::Zero());
    ControlVectorArray<control_dim> ff_traj(K_, ct::core::ControlVector<control_dim>::Zero());

    NLOptConSolver_t::Policy_t init_controller(state_ref_traj, ff_traj, fb_traj, settings_.dt);
    solver_ptr_->setInitialGuess(init_controller);
    solver_ptr_->solve();
    std::cout << "In constructor of the MPC class" << std::endl;;
    std::cout << "Initial solving OK" << std::endl;

    StateFeedbackController<state_dim, control_dim> init_solution = solver_ptr_->getSolution();


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
           std::shared_ptr<StateFeedbackController<state_dim, control_dim>> new_policy_ptr)
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

    bool success = mpc_ptr_->finishIteration(
        estimated_state, t, *new_policy_ptr, ts_new_policy);

    return success;
  }
};


double constrainAngle(double x){
    x = fmod(x + M_PI, 2*M_PI);
    if (x < 0)
        x += 2*M_PI;
    return x - M_PI;
}


int main(int argc, char **argv)
{
  using Scalar = ct::core::ADCGScalar;
  const Scalar M(0.5);
  const Scalar m(0.2);
  const Scalar b(0.1);
  const Scalar l(0.3);
  const Scalar I(0.006);


  ros::init(argc, argv, "mpc_demo_node");
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

  // std::shared_ptr<InvpendSystem<Scalar>> invpend_sys_ptr(
      // new InvpendSystem<Scalar>(M, m, b, l, I));

  std::shared_ptr<InvpendSystem<double>> invpend_sys_double_ptr(
      new InvpendSystem<double>(0.5, 0.2, 0.1, 0.3, 0.006));


  ControlVector<control_dim> u_max, u_min;
  u_max << 5.0;
  u_min << -5.0;
  std::string load_path = "/home/baozhe/FastLab/report_ws/src/control_demo/parameters";


  InvpendMPC mpc(load_path, u_max, u_min, invpend_sys_double_ptr);

  size_t window_num = mpc.K_ + 1;
  StateVectorArray<state_dim> ref;
  for (size_t i = 0; i < window_num; i++) {
    ref.push_back(StateVector<state_dim>::Zero());
  }
  auto new_policy_ptr = std::make_shared<StateFeedbackController<state_dim, control_dim>>();
  
  ros::Rate freq(100.0);
  auto start = ros::Time::now().toSec();
  while (ros::ok() && mpc.run(CURRENT_EST, ref, new_policy_ptr)) {
    // ROS_INFO("EST: angle: %.2f", CURRENT_EST(1));

    auto t = ros::Time::now().toSec() - start;

    std_msgs::Float64 msg;
    auto command = new_policy_ptr->uff()[0];
    msg.data = command(0);
    effort_pub.publish(msg);

    freq.sleep();
    ros::spinOnce();
  }


  return 0;

}