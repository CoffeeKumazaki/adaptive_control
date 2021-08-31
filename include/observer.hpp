#pragma once

template<class state, unsigned char delimiter>
struct Observer{
  std::ofstream fout;
  Observer(const std::string& FileName) :fout(FileName){};

  void operator()(const state& x, double t){
    fout << t;
    for (size_t i = 0; i < x.size(); i++) {
      fout << delimiter << x[i]; 
    }
    fout << std::endl;
  }
};