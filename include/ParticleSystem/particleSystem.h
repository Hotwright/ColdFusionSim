#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H

const uint64_t NULL_INDEX = uint64_t(-1);

#include <vector>
#include <time.h>
#include <math.h>
#include <random>
#include <iostream>
#include <thread>

std::default_random_engine generator;
std::uniform_real_distribution<double> U(0.0,1.0);
std::normal_distribution<double> normal(0.0,1.0);

class ParticleSystem{

  friend class ParticleSystemRenderer;
  friend class Trajectory;

public:

  ParticleSystem(
    uint64_t N,
    uint64_t capacity,
    double dt = 1.0/300.0,
    double density = 0.25,
    double Lx = 0.5, double Ly = 1.0, double Lz = 0.5,
    uint64_t seed = clock()
  )
  // radius chosen so N spheres of this radius occupy `density` of the box volume.
  //  N is the Pd/Ni cap (matches the Particles slider); capacity is the total
  //  storage reserved, since proteium (H) particles are added on top of N.
  : nParticles(N), maxCapacity(capacity), radius(std::cbrt(3.0*density*Lx*Ly*Lz/(4.0*M_PI*N))),drag(0),
    mass(1.0),dt(dt),collisionTime(10*dt),
    shakerPeriod(1.0),shakerAmplitude(radius),shakerTime(0.0),
    Lx(Lx), Ly(Ly), Lz(Lz),
    pdAttraction(0.0),attractionRangeMultiplier(2.5),
    hRadiusRatio(1.0),hRadius(radius),hMass(mass),nH(0),
    // formal charges (simulation units, not SI): Pd2+/Ni2+ vs hydride H-,
    //  so proteium is drawn toward the metal and metal ions repel each other
    pdCharge(2.0),niCharge(2.0),hCharge(-1.0),coulombStrength(3.0)
  {

    floatState = new float [capacity*4];

    generator.seed(seed);
    Nc = std::ceil(1.0/(4.0*radius));
    deltax = Lx / Nc;
    deltay = Ly / Nc;
    deltaz = Lz / Nc;

    shakerDisplacement = shakerAmplitude/2.0;
    massRatio = 1.0;
    radiusRatio = 2.0;

    for (int c = 0; c < Nc*Nc*Nc; c++){
      cells.push_back(NULL_INDEX);
    }

    for (int i = 0; i < capacity; i++){
      list.push_back(NULL_INDEX);
    }

    for (int i = 0; i < N; i++){
      double x = U(generator)*(Lx-2*radius)+radius;
      double y = U(generator)*(Ly-2*radius)+radius;
      double z = U(generator)*(Lz-2*radius)+radius;

      double r = i%2 == 0 ? radius : radius / radiusRatio;
      double m = i%2 == 0 ? mass : mass / massRatio;
      double q = i%2 == 0 ? pdCharge : niCharge;

      addParticle(x,y,z,r,m,q);
      uint64_t c = hash(i);
      if (cells[c] == NULL_INDEX){
        cells[c] = i;
      }
      else{
        insert(cells[c],uint64_t(i));
      }
    }

    setCoeffientOfRestitution(0.95);
  }

  void applyForce(double fx, double fy, double fz);

  void step();

  // grows/shrinks the Pd/Ni population, capped at N; inserts/removes
  //  before the trailing proteium (H) block so H particles are untouched
  void addParticle();
  void removePdNiParticle(){removeParticle(size()-1-nH);}

  void removeParticle(){removeParticle(size()-1);}

  uint64_t size(){return uint64_t(std::floor(state.size() / 3));}

  // for judging separation, 1 means perfect separation, 0 means random
  //   particle placement, works 'better' the more particles there are
  double orderParameter();

  // parameter setters

  void randomise(double propBig);

  void setTimeStep(double dt){ if(this->dt!=dt) {newTimeStepStates(this->dt,dt);} this->dt = dt; dtdt = dt*dt; }

  void setShakerPeriod(double p){
    if (p != shakerPeriod) {shakerTime = shakerTime*p/shakerPeriod;}
    shakerPeriod = p;
  }

  void setMassRatio(double mr){
    if (massRatio != mr){massRatio = mr; randomise(0.5);}
  }

  void setRadiusRatio(double rr){
    if (radiusRatio != rr){radiusRatio = rr; randomise(0.5);}
  }

  // attractive force between big (Pd) particles, active from contact
  //  out to attractionRangeMultiplier times the contact distance
  void setPdAttraction(double a){pdAttraction = a;}
  double getPdAttraction(){return pdAttraction;}

  // fixes the proteium (H) particle size, given as the ratio of the
  //  Pd (base) radius to the H radius; mass scales with volume so H
  //  particles stay light, matching hydrogen's low atomic mass
  void setHRadiusRatio(double r){
    hRadiusRatio = r;
    hRadius = radius / hRadiusRatio;
    hMass = mass * (hRadius/radius) * (hRadius/radius) * (hRadius/radius);
  }

  // adds/removes proteium (H) particles, on top of the Pd/Ni population,
  //  to reach the target count; no-op if already there
  void setHCount(uint64_t target);
  uint64_t getHCount(){return nH;}

  void setShakerAmplitude(double a){
    if (shakerAmplitude != a*radius){
      shakerAmplitude = a*radius; shakerDisplacement = shakerAmplitude;
    }
  }

  // depends on collision time (which should be >> timestep, 10x seems to work)
  //  and the masses of two colliding particles.
  void setCoeffientOfRestitution(double c);

  // helper for setting restitution
  double reducedMass(float m1, float m2){
      return 1.0 / (1.0/m1 + 1.0/m2);
  }

  // helper for setting restitution
  double damping(float m1, float m2){
      double meff = reducedMass(m1,m2);
      return 2*meff*(-std::log(coefficientOfRestitution)/collisionTime);
  }

  // helper for setting restitution
  double restoration(float m1, float m2){
      double meff = reducedMass(m1,m2);
      return meff/(collisionTime*collisionTime)*(std::log(coefficientOfRestitution)+M_PI*M_PI);
  }

  // parameter getters

  double getshakerPeriod(){return shakerPeriod;}
  double getShakerAmplitude(){return shakerAmplitude/radius;}
  double getLx(){return Lx;}
  double getLy(){return Ly;}
  double getLz(){return Lz;}

  // scenarios

  void one();

  // scatters every particle to a random position in the box, zeroing
  //  velocity; species/radius/mass are untouched
  void randomisePositions();

  ~ParticleSystem(){
    free(floatState);
  }

private:

  std::vector<double> state;
  std::vector<double> lastState;

  std::vector<double> parameters;
  std::vector<double> charge;

  std::vector<double> forces;
  std::vector<double> velocities;

  std::vector<uint64_t> cells;
  std::vector<uint64_t> list;

  double Lx;
  double Ly;
  double Lz;

  uint64_t Nc;
  double deltax;
  double deltay;
  double deltaz;

  uint64_t nParticles;
  uint64_t maxCapacity;

  double coefficientOfRestitution;
  double collisionTime;

  double massRatio;
  double radiusRatio;

  double radius;
  double drag;
  double mass;
  double dt;
  double dtdt;

  double shakerDisplacement;
  double shakerPeriod;
  double shakerAmplitude;
  double shakerTime;

  double pdAttraction;
  double attractionRangeMultiplier;

  double hRadiusRatio;
  double hRadius;
  double hMass;
  uint64_t nH;

  double pdCharge;
  double niCharge;
  double hCharge;
  double coulombStrength;

  float * floatState;

  // appends at the true end of the arrays (used for H particles, which
  //  are always kept as the trailing block)
  void addParticle(double x, double y, double z, double r, double m, double q);
  // inserts a Pd/Ni particle just before the trailing H block
  void insertParticle(uint64_t idx, double x, double y, double z, double r, double m, double q);
  void removeParticle(uint64_t i);

  // Cell Linked List Collisions detection
  void resetLists();
  void insert(uint64_t next, uint64_t particle);
  void populateLists();
  void handleCollision(uint64_t i, uint64_t j);
  void cellCollisions(
    int64_t a1, int64_t b1, int64_t c1,
    int64_t a2, int64_t b2, int64_t c2
  );

  uint64_t hash(float x, float y, float z){
    return uint64_t(floor(x/deltax))*Nc*Nc + uint64_t(floor(y/deltay))*Nc + uint64_t(floor(z/deltaz));
  }

  uint64_t hash(uint64_t particle){
    uint64_t h = uint64_t(floor(state[particle*3]/deltax))*Nc*Nc
      + uint64_t(floor(state[particle*3+1]/deltay))*Nc
      + uint64_t(floor(state[particle*3+2]/deltaz));
    // a particle outside the box would normally crash the program
    //  (segfault), if this happens there is no logic to get it back in currently
    //  these cases happen due to numerical instability, so a smaller time step
    //  will help! - pretty rare with defaults
    if (h < 0 || h > Nc*Nc*Nc){return 0;}
    return h;
  }

  // integrator requires current and last positions
  //  changing the timestep can introduce instability if
  //  not accounted for
  void newTimeStepStates(double oldDt, double newDt);
};

#endif
