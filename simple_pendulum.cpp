/**
 * @file simple_pendulum.cpp
 * @brief Main file for numerical simulation of "wall hop" modeled as simple
 * pendulum.
 * @author Takanishi Labratory, Environmental Monitoring team, Jiei Suzuki.
 * @date 2018/11/17
 */

#include <cmath>
#include <iostream>
#include <vector>

#include "matplotlib-cpp-master/matplotlibcpp.h" //グラフの描画

namespace plt = matplotlibcpp;

struct Point2D {
  double x;
  double y;
};

struct condition_param { //初期条件とかパラメータ
  double M;              //[kg]
  double L;              //[m]
  double theta;          //[deg]
  double v_0;            //[m/s]
  double g;              //[m/s^2]
};

class simple_pendulum {
public:
  double dT;
  double
      time_limit; //タイムリミットになったら終了？or速度がマイナスになったら終了？

  std::vector<double> time;
  std::vector<Point2D> abs_mass_point; //質点の座標

  double x_2dot, x_dot;
  double phi_dot, phi;
  condition_param pc;

  Point2D point_tmp;

  simple_pendulum(condition_param input_pc) {
    pc = input_pc;
    dT = 0.001;
    time_limit = 5.0;

    x_2dot = 0;
    x_dot = pc.v_0;
    phi_dot = 0;
    phi = 0;

    time.push_back(0.0);

    point_tmp.x = (-1) * pc.L * cos(pc.theta / 180 * M_PI);
    point_tmp.y = (-1) * pc.L * sin(pc.theta / 180 * M_PI);
    abs_mass_point.push_back(point_tmp);
  }

  void renew_state() { //考えるのめんどいからオイラー法
    x_2dot = pc.g * cos(pc.theta / 180 * M_PI - phi);
    x_dot -= dT * x_2dot;
    phi_dot = x_dot / pc.L;
    phi += dT * phi_dot;
    time.push_back(time.back() + dT);
    point_tmp.x = (-1) * pc.L * cos(pc.theta / 180 * M_PI - phi);
    point_tmp.y = (-1) * pc.L * sin(pc.theta / 180 * M_PI - phi);
    abs_mass_point.push_back(point_tmp);
    std::cout << "time=" << time.back() << "\tx=" << point_tmp.x
              << "\ty=" << point_tmp.y << '\n';
    // std::cout << "x_dot=" << x_dot << '\n';
  }
};

int main(int argc, char const *argv[]) {
  //条件入力
  condition_param pc = {3.0, 1.0, 80, 1.0, 9.8};
  simple_pendulum sp(pc);

  while ((sp.time.back() < sp.time_limit) && (sp.x_dot > 0)) {
    sp.renew_state();
  }

  /*std::cout << "H = " << sp.pc.L * sin(sp.phi) << '\n';
  std::vector<double> x_tmp(1), y_tmp(1);
  for (unsigned int i = 0; (i * 10) < sp.abs_mass_point.size(); i++) {
    x_tmp[0] = sp.abs_mass_point[i * 10].x;
    x_tmp[0] = sp.abs_mass_point[i * 10].y;
    plt::plot(x_tmp, y_tmp, "xr");
    plt::pause(sp.dT);
  }*/
  return 0;
}
