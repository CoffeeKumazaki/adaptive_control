#include <stdafx.hpp>
#include <system/LTI_system.hpp>

using namespace std;

Eigen::MatrixXd SolveRiccati(
  const Eigen::MatrixXd &A,
  const Eigen::MatrixXd &B,
  const Eigen::MatrixXd &Q,
  const Eigen::MatrixXd &R) 
{
  const uint dim_x = A.rows();
  const uint dim_u = B.cols();

  // set Hamilton matrix
  Eigen::MatrixXd Ham = Eigen::MatrixXd::Zero(2 * dim_x, 2 * dim_x);
  Ham << A, -B * R.inverse() * B.transpose(), -Q, -A.transpose();

  // calc eigenvalues and eigenvectors
  Eigen::EigenSolver<Eigen::MatrixXd> Eigs(Ham);

  // check eigen values
  // std::cout << "eigen values：\n" << Eigs.eigenvalues() << std::endl;
  // std::cout << "eigen vectors：\n" << Eigs.eigenvectors() << std::endl;

  // extract stable eigenvectors into 'eigvec'
  Eigen::MatrixXcd eigvec = Eigen::MatrixXcd::Zero(2 * dim_x, dim_x);
  int j = 0;
  for (int i = 0; i < 2 * dim_x; ++i) {
    if (Eigs.eigenvalues()[i].real() < 0.) {
      eigvec.col(j) = Eigs.eigenvectors().block(0, i, 2 * dim_x, 1);
      ++j;
    }
  }

  // calc P with stable eigen vector matrix
  Eigen::MatrixXcd Vs_1, Vs_2;
  Vs_1 = eigvec.block(0, 0, dim_x, dim_x);
  Vs_2 = eigvec.block(dim_x, 0, dim_x, dim_x);
  Eigen::MatrixXd P = (Vs_2 * Vs_1.inverse()).real();

  return P;
}

class MRAC_SISO {
  static const int PARAM_NUM = 2+1;
public:
  using parameter_t = Eigen::Vector<double, PARAM_NUM>;
  using state_t  = LTISystem<1, 2, 2>::state_t;
  using input_t  = LTISystem<1, 2, 2>::input_t;
  using output_t = LTISystem<1, 2, 2>::output_t;

public:
  MRAC_SISO(state_t init) 
  : model("model", init)
  {
    gain_ky << 1.0, 0.0, 
               0.0, 1.0;
    gain_kr = 1.0;

    model.a << 0.0,  1.0, 
              -1.0, -1.0;
    model.b << 0.0, 1.0;
    model.c << 1.0, 0.0, 
               0.0, 1.0;
    model.d << 0.0, 0.0;

    parameter << 0.01, 0.01, 0.01;

    P = SolveRiccati(model.a, Eigen::MatrixXd::Zero(2,1), Eigen::MatrixXd::Identity(2,2), Eigen::MatrixXd::Identity(1,1));
  }
  ~MRAC_SISO() {}

public:
  double decision(const output_t& y, const input_t& r, double dt) {

    model.input = r;
    model.update(dt);

    ym = model.observation(); 

    output_t e = ym - y;

    kd << gain_ky * y * e.transpose() * P.col(1),
          gain_kr * r * e.transpose() * P.col(1);

    boost::numeric::odeint::integrate_const(
      stepper, *this, parameter, 0.0, dt, dt*0.1
    );

    return parameter[0]*y[0] + parameter[1]*y[1] + parameter[2]*r[0];
  }

  void operator()(const parameter_t&x, parameter_t& dx, double dt) {

    dx = kd;
  }

  parameter_t parameter;
  output_t ym;

private:
  LTISystem<1, 2, 2> model;
  Eigen::Matrix<double, 2, 2> P;
  Eigen::Matrix<double, 2, 2> gain_ky;
  double gain_kr;
  parameter_t kd; 
  boost::numeric::odeint::runge_kutta4<parameter_t> stepper;
};

int main(int argc, char const *argv[])
{
  LTISystem<1, 2, 2>::state_t init;
  init << 0.0, 0.0;

  LTISystem<1, 2, 2> plant("plant", init);
  plant.a << 0.0,  1.0, 
            -4.0, -2.0;
  plant.b << 0.0, 4.0;
  plant.c << 1.0, 0.0, 
             0.0, 1.0;
  plant.d << 0.0, 0.0;
  
  MRAC_SISO mrac(init);

  double t = 0;
  double dt = 0.1;
  double simTime = 100.0;

  while(t < simTime) {
    
    MRAC_SISO::output_t y = plant.observation();
    MRAC_SISO::input_t input; input << sin(2.0*t);
    double mrac_input = mrac.decision(y, input, dt);

    plant.input << mrac_input;

    plant.update(dt);

    cout << t << " " << y[0] << " " << y[1] << " "
              << mrac.ym[0] << " " << mrac.ym[1] << " "
              << mrac_input << " " 
              << mrac.parameter[0] << " " << mrac.parameter[1] << " " << mrac.parameter[2] << endl;

    t += dt;
  }

  return 0;
}