#include <stdafx.hpp>
#include <filesystem>
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

public:
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
    if (t >= (n-m)) {
      for (size_t i = 0; i < n && i <= t-(n-m); i++) {
        phit[i] = state_history[t-(n-m)-i][X];
      }
      for (int i = 1; i < n; i++) {
        phit[n+i] = input_history[t-(n-m)-i];
      }
    }

    // calc error.
    double error = 0.0;

    if (parameters.size() > 0) {
      parameter prevParam = parameters.back();
      error = yt[X] - prevParam.transpose() * phit;
    }
    std::cout << error << std::endl;

    // update parameter.
    parameter pt = parameter::Random();
    if (parameters.size() > 0) {
      pt = parameters.back() + gain * phit / ( c + phit.transpose()*phit) * error;
    }

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

void run_simulation(std::string name, double k, double cs, double a, double c) {
  
  state x0;
  x0 << 0.0, 0.0;

  System ref("reference", x0);
  System plant("plant", x0);
  Observer<state, ' '> refObs(name + "_mrac_ref.dat");
  Observer<state, ' '> pltObs(name + "_mrac_plt.dat");
  Observer<MRAC<2,0>::parameter, ' '> prmObs(name + "_mrac_prm.dat");

  plant.a << 0.0, 1.0, -k, -cs;

  double t = 0;
  double dt = 0.1;
  double simTime = 100.0;

  MRAC<2,0> mrac(a, c);
  mrac.ref = &ref;

  while(t < simTime) {
    
    double input = step(t);

    ref.input << input;

    double mrac_input = mrac.decision(plant.state, input, dt);
    plant.input << mrac_input;

    ref.update(dt);
    plant.update(dt);

    t += dt;
    refObs(ref.state, t);
    pltObs(plant.state, t);
    prmObs(mrac.parameters.back(), t);  
  }

}


int main(int argc, char const *argv[])
{
  if (argc <= 1) {
    cout << "no file" << endl;
    return 1;
  }
  std::filesystem::path filename = argv[1];

  string dirname = filename.filename();
  std::filesystem::create_directory(dirname);

  ifstream ifs(filename);
  if (!ifs) {
    cout << "file not found" << endl;
    return 1;
  }

  while (!ifs.eof()) {
    string name;
    double k, cs, a, c;

    ifs >> name >> k >> cs >> a >> c;
    cout << "run " + name << endl;
    run_simulation(dirname + "/"  + name, k, cs, a, c);
  }

  return 0;
}