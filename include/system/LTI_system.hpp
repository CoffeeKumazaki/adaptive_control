#pragma once

template <size_t input_dim, size_t output_dim, size_t state_dim>
class LTISystem {

public:
  using state_t = Eigen::Vector<double, state_dim>;
  using input_t = Eigen::Vector<double, input_dim>;
  using output_t = Eigen::Vector<double, output_dim>;

public:
  LTISystem(std::string name, const state_t& init) 
  : t(0)
  {
    state = init;
    input = input_t::Zero();
  }
  ~LTISystem(){}

  void update(double dt) {

    boost::numeric::odeint::integrate_const(
      stepper, *this, state, t, t+dt, dt*0.1
    );
  }

  output_t observation() {

    return c * state + d * input;
  }

  void operator()(const state_t&x, state_t& dx, double dt) {

    dx = a * x + b * input;
  }

public:
  double t;
  boost::numeric::odeint::runge_kutta4<state_t> stepper;

  state_t state;
  input_t input;
  Eigen::Matrix<double, state_dim,  state_dim> a;
  Eigen::Matrix<double, state_dim,  input_dim> b;
  Eigen::Matrix<double, output_dim, state_dim> c;
  Eigen::Matrix<double, output_dim, input_dim> d;
};

template <size_t input_dim, size_t output_dim, size_t state_dim>
std::ostream& operator<<(std::ostream& os, const LTISystem<input_dim, output_dim, state_dim>& sys) {

  os << "# a --" << std::endl << sys.a << std::endl << std::endl;
  os << "# b --" << std::endl << sys.b << std::endl << std::endl;
  os << "# c --" << std::endl << sys.c << std::endl << std::endl;
  os << "# d --" << std::endl << sys.d << std::endl << std::endl;
  os << "# state --" << std::endl << sys.state << std::endl << std::endl;
  os << "# input --" << std::endl << sys.input << std::endl << std::endl;

  return os;
}
