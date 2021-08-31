#include <stdafx.hpp>
#include <observer.hpp>
#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>

using namespace std;

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
  System(std::string name, const state& init)
  {
    state = init;

    a << 0.0, 1.0, -0.2, -0.5;
    b << 0.0, 1.0;
    c << 1.0, 0.0, 0.0, 1.0;
    input << 0.0;
  }

  void update(double dt) {

    boost::numeric::odeint::integrate_const(
      stepper, *this, state, t, t+dt, dt*0.1
    );
  }

  void operator()(const state&x, state& dx, double dt) {

    dx = a * x + b * input;
  }

public:
  double t;
  boost::numeric::odeint::runge_kutta4<state> stepper;

  state state;
  Eigen::Matrix<double, STATE_NUM, STATE_NUM> a;
  Eigen::Matrix<double, STATE_NUM, INPUT_NUM> b;
  Eigen::Matrix<double, STATE_NUM, STATE_NUM> c;
  input input;

};

// model reference adaptive control
#define N (2)
#define M (0)
#define D (N-M)

template <size_t n, size_t m>
class MRAC {

  using parameter = Eigen::Vector<double, (n-m)*2>;
  using phi = Eigen::Vector<double, (n-m)*2>;

public:
  MRAC(double gain_, double c_) : gain(gain_), c(c_) {
    state_history.clear();
    input_history.clear();
  }

  double decision(const state& yt, double r, double dt) {

    state_history.push_back(yt);

    phi phit = phi::Zero();
    int t = state_history.size() - 1;
    cout << "t: " << t  << " yt: " << yt << endl;
    if (t >= (n-m)) {
      for (size_t i = 0; i < n && i <= t-(n-m); i++) {
        phit[i] = state_history[t-(n-m)-i][X];
      }
      for (int i = 1; i < n; i++) {
        phit[n+i] = input_history[t-(n-m)-i];
      }
    }
    cout << "phi(t-d)" << endl;
    cout << phit << endl;

    // calc error.
    double error = 0.0;

    if (parameters.size() > 0) {
      parameter prevParam = parameters.back();
      error = yt[X] - prevParam.transpose() * phit;
    }
    cout << "error " << error << endl;

    // update parameter.
    parameter pt = parameter::Random();
    if (parameters.size() > 0) {
      pt = parameters.back() + gain * phit / ( c + phit.transpose()*phit) * error;
    }
    cout << "pt" << endl;
    cout << pt << endl;

    parameters.push_back(pt);

    // calc input.
    double u = r;
    if (t >= (n-m)) {
      u = 0.0;
      for (size_t i = 0; i < n; i++) {
        u += -pt[i]*state_history[t-i][X];
      }
      for (size_t i = 1; i < n; i++) {
        u += -pt[n+i]*input_history[t-i];
      }

      System rref = *ref;
      rref.input << r;
      rref.update(dt*(n-m));
      double ymd = rref.state[X];
      u += ymd;
      
      u /= pt[n];
    }
    cout << "u: " << u << " r: " << r << endl;

    input_history.push_back(u);

    return u;
  }

  double a(double t) { return 1.0; }
public:
  System* ref;
  double c;
  double gain;
  std::vector<parameter> parameters;
  std::vector<state>  state_history;
  std::vector<double> input_history;
};


int main(int argc, char const *argv[])
{
  state x0;
  x0 << 0.0, 0.0;

  System ref("reference", x0);
  System plant("plant", x0);
  Observer<state, ' '> refObs("mrac_ref.dat");
  Observer<state, ' '> pltObs("mrac_plt.dat");

  plant.a << 0.0, 1.0, -0.2, -0.2;

  double t = 0;
  double dt = 0.1;
  double simTime = 200.0;

  MRAC<2,0> mrac(1.2, 0.1);
  mrac.ref = &ref;

  while(t < simTime) {
    std::cout << "t:  " << t << std::endl;
    
    double input = sinwave(t);

    ref.input << input;

    std::cout << " decision ..." << std::endl;
    double mrac_input = mrac.decision(plant.state, input, dt);
    std::cout << "mrac_input " << mrac_input << std::endl;
    plant.input << mrac_input;

    ref.update(dt);
    plant.update(dt);
    std::cout << "plant " << plant.state[X] << std::endl;

    t += dt;
    refObs(ref.state, t);
    pltObs(plant.state, t);
    std::cout << std::endl;
  }

  return 0;
}
