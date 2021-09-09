#include <stdafx.hpp>
#include <system/LTI_system.hpp>

using namespace std;

class MRAC_SISO {
  static const int PARAM_NUM = 2;
  using parameter_t = Eigen::Vector<double, PARAM_NUM>;
  using state_t  = LTISystem<1, 1, 1>::state_t;
  using input_t  = LTISystem<1, 1, 1>::input_t;
  using output_t = LTISystem<1, 1, 1>::output_t;

public:
  MRAC_SISO(state_t init) 
  : model("model", init)
  {
    gain_a = 1.0;
    gain_b = 1.0;

    model.a << -1.0;
    model.b << 1.0;
    model.c << 1.0;
    model.d << 0.0;

    parameter << 0.1, 1;
  }
  ~MRAC_SISO() {}

public:
  double decision(double y, double r, double dt) {

    model.input << r;
    model.update(dt);

    ym = model.observation()[0]; 

    double e = ym - y;
    double ky = calc_ky();
    double kr = calc_kr();
    double u = ky * y + kr * r;

    pd << -gain_a * y * e, 
          -gain_b * u * e;

    boost::numeric::odeint::integrate_const(
      stepper, *this, parameter, 0.0, dt, dt*0.1
    );

    return calc_ky()*y + calc_kr()*r;
  }

  double calc_ky() {
    return (model.a[0] - parameter[0])/parameter[1];
  }

  double calc_kr() {
    return model.b[0]/parameter[1];
  }


  void operator()(const parameter_t&x, parameter_t& dx, double dt) {

    dx = pd;
  }

  parameter_t parameter;
  double ym;

private:
  LTISystem<1, 1, 1> model;
  double gain_a, gain_b;
  parameter_t pd; 
  boost::numeric::odeint::runge_kutta4<parameter_t> stepper;
};

int main(int argc, char const *argv[])
{
  LTISystem<1, 1, 1>::state_t init;
  init << 0.0;

  LTISystem<1, 1, 1> plant("plant", init);
  plant.a << 0.5;
  plant.b << 1.0;
  plant.c << 1.0;
  plant.d << 0.0;
  
  MRAC_SISO mrac(init);

  double t = 0;
  double dt = 0.1;
  double simTime = 100.0;

  while(t < simTime) {
    
    double y = plant.observation()[0];
    double input = sin(t);
    double mrac_input = mrac.decision(y, input, dt);

    plant.input << mrac_input;

    plant.update(dt);

    cout << t << " " << y << " " << mrac.ym << " " << mrac_input << " " << mrac.parameter[0] << " " << mrac.parameter[1] << endl;

    t += dt;
  }

  return 0;
}
