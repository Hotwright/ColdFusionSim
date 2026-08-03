#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <ParticleSystem/particleSystem.h>

#include <fstream>
#include <utils.h>
#include <ctime>
#include <thread>

struct State {
  State (double x, double y, double z, double r)
  : x(x), y(y), z(z), radius(r)
  {}
  double x;
  double y;
  double z;
  double radius;
};

struct Parameters {
  Parameters(double m, double r, double p, double a, double re, double n, double o, double s)
  : massRatio(m), radiusRatio(r), period(p), amplitude(a), restitution(re), nParticles(n),order(o),speed(s)
  {}
  double massRatio;
  double radiusRatio;
  double period;
  double amplitude;
  double restitution;
  double nParticles;
  double order;
  double speed;
};

class Trajectory {
public:

const int CACHE_SIZE_BYTES = 1e6;

  Trajectory(){
    newFile();
  }

  void takeReading(ParticleSystem & p);

  void setSpeed(double s){speed = s;}

  void save();

  void clear(){trajectory.clear(); shakerPosition.clear(); parameters.clear(); size = 0;}

  void newFile(){
    clear();
    ioMode = std::ios_base::trunc;
    std::time_t t = std::time(NULL);
    std::tm tm = *std::localtime(&t);
    std::stringstream ss;
    ss << "trajectory-";
    ss << std::put_time(&tm, "%c %Z");
    file = ss.str();
    size = 0;
  }

  std::string fileName(){return file;}

private:

  // timestep implicitly defined by use of takeReading()
  std::vector<std::vector<State>> trajectory;
  std::vector<double> shakerPosition;
  std::vector<Parameters> parameters;
  std::ios_base::openmode ioMode;
  std::string file;
  double speed;

  uint64_t size;

  void threadedSave(std::vector<std::vector<State>> trajectory);
};

#endif
