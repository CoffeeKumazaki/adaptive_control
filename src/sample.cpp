#include <stdafx.hpp>
#include <observer.hpp>
#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>

enum STATE {
  X = 0,
  VX,
  STATE_NUM,
};

enum INPUT {
  F = 0,
  INPUT_NUM,
};

using state = Eigen::Vector<double, STATE_NUM>;
using input = Eigen::Vector<double, INPUT_NUM>;

double step(double t) {

  return (t<0.0) ? 0.0 : 1.0;
}

#define omega (0.05)
double sinwave(double t) {

  return sin(omega*t);
}

typedef double (*InputFunc)(double); 

class System {

public:
  System(const state& init) {
    state = init;

    a << 0.0, 1.0, -0.2, -0.5;
    b << 0.0, 1.0;
    c << 1.0, 0.0, 0.0, 1.0;
    input << 0.0;
  }

  void operator()(const state&x, state& dx, double dt) {

    dx = a * x + b * input;
  }

public:
  state state;
  Eigen::Matrix<double, STATE_NUM, STATE_NUM> a;
  Eigen::Matrix<double, STATE_NUM, INPUT_NUM> b;
  Eigen::Matrix<double, STATE_NUM, STATE_NUM> c;
  input input;
};

int main(int argc, char const *argv[])
{
  state x0;
  x0 << 0.0, 0.0;

  System ref(x0);
  System plant(x0);
  plant.a << 0.0, 1.0, -1.2, -1.2;
  Observer<state, ' '> observer_ref("ref.dat");
  Observer<state, ' '> observer_plt("plant.dat");

  boost::numeric::odeint::runge_kutta4<state> stepper;

  double t = 0;
  double dt = 0.1;
  double simTime = 100.0;

  while(t < simTime) {

    observer_ref(ref.state, t);
    observer_plt(plant.state, t);

    double input = sinwave(t);

    ref.input << input;
    plant.input << input;
    boost::numeric::odeint::integrate_const(
      stepper, ref, ref.state, t, t+dt, dt*0.1
    );
    boost::numeric::odeint::integrate_const(
      stepper, plant, plant.state, t, t+dt, dt*0.1
    );
    t += dt;
  }

  return 0;
}
