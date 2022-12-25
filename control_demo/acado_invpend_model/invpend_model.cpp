#include <memory>
#include <acado_optimal_control.hpp>
#include <acado_code_generation.hpp>
#include <acado_gnuplot.hpp>

int main() {
  USING_NAMESPACE_ACADO

  DifferentialState     x, x_dot, theta, theta_dot;
  Control               u;
  DifferentialEquation  f;
  Function              h, hN;

  const double M        = 0.5;       // cart mass
  const double m        = 0.2;       // pole mass
  const double b        = 0.1;       // friction parameter
  const double l        = 0.3;       // half len of pole
  const double I        = 0.006;     // mass inertia moment
  const double G        = 9.8;     
  const double t_start  = 0.0;
  const double t_end    = 2.0;
  const double dt       = 0.01;
  const int N           = round(t_end / dt);
  const double u_max    = 2.0;
  const double u_min    = -2.0;

  f << dot(x)         == x_dot;
  f << dot(x_dot)     == 1 / ( M + m - (m * m * l * l * cos(theta) * cos(theta)) / (I + m * l * l) ) * 
        ( u - b * x_dot + m * l * sin(theta) * theta_dot * theta_dot - m * l * cos(theta) * m * G * l * sin(theta) / (I + m * l * l) );
  f << dot(theta)     == theta_dot;
  f << dot(theta_dot) == 1 / (I + m * l * l) * 
        (m * G * l * sin(theta) - m * l * (1 / ( M + m - (m * m * l * l * cos(theta) * cos(theta)) / (I + m * l * l) ) * 
        ( u - b * x_dot + m * l * sin(theta) * theta_dot * theta_dot - m * l * cos(theta) * m * G * l * sin(theta) / (I + m * l * l) )) * cos(theta));


  h << x << x_dot << theta << theta_dot << u;
  hN << x << x_dot << theta << theta_dot;

  OCP ocp(t_start, t_end, N);

  // For code generation, references are set during run time.
  BMatrix Q_sparse(h.getDim(), h.getDim());
  Q_sparse.setIdentity();
  BMatrix QN_sparse(hN.getDim(), hN.getDim());
  QN_sparse.setIdentity();
  ocp.minimizeLSQ( Q_sparse, h);
  ocp.minimizeLSQEndTerm( QN_sparse, hN );

  ocp.subjectTo(f);
  ocp.subjectTo(u_min <= u <= u_max);

  // For code generation, we can set some properties.
  // The main reason for a setting is given as comment.
  OCPexport mpc(ocp);

  mpc.set(HESSIAN_APPROXIMATION,  GAUSS_NEWTON);        // is robust, stable
  mpc.set(DISCRETIZATION_TYPE,    MULTIPLE_SHOOTING);   // good convergence
  mpc.set(SPARSE_QP_SOLUTION,     FULL_CONDENSING_N2);  // due to qpOASES
  mpc.set(INTEGRATOR_TYPE,        INT_IRK_GL4);         // accurate
  mpc.set(NUM_INTEGRATOR_STEPS,   N);
  mpc.set(QP_SOLVER,              QP_QPOASES);          // free, source code
  mpc.set(HOTSTART_QP,            YES);
  mpc.set(CG_USE_OPENMP,                    YES);       // paralellization
  mpc.set(CG_HARDCODE_CONSTRAINT_VALUES,    NO);        // set on runtime
  mpc.set(CG_USE_VARIABLE_WEIGHTING_MATRIX, YES);       // time-varying costs
  mpc.set( USE_SINGLE_PRECISION,        YES);           // Single precision

  // Do not generate tests, makes or matlab-related interfaces.
  mpc.set( GENERATE_TEST_FILE,          NO);
  mpc.set( GENERATE_MAKE_FILE,          NO);
  mpc.set( GENERATE_MATLAB_INTERFACE,   NO);
  mpc.set( GENERATE_SIMULINK_INTERFACE, NO);

  // Finally, export everything.
  if(mpc.exportCode("invpend_mpc_codegen") != SUCCESSFUL_RETURN)
    exit( EXIT_FAILURE );
  mpc.printDimensionsQP( );

  return EXIT_SUCCESS;
}