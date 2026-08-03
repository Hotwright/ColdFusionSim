#include <ParticleSystem/particleSystem.h>
#include <time.h>

void ParticleSystem::resetLists(){
  for (int i = 0; i < Nc*Nc*Nc; i++){
    cells[i] = NULL_INDEX;
  }
  // list is indexed by particle slot, sized to the full capacity since
  //  proteium (H) particles can push the active count past nParticles
  for (int i = 0; i < list.size(); i++){
    list[i] = NULL_INDEX;
  }
}

void ParticleSystem::insert(uint64_t next, uint64_t particle){
  uint64_t i = next;
  while (list[i] != NULL_INDEX){
    i = list[i];
  }
  list[i] = particle;
}

void ParticleSystem::populateLists(
){
  for (int i = 0; i < size(); i++){
    uint64_t c = hash(i);
    if (cells[c] == NULL_INDEX){
      cells[c] = uint64_t(i);
    }
    else{
      insert(cells[c],uint64_t(i));
    }
  }
}

void ParticleSystem::handleCollision(uint64_t i, uint64_t j){
  if (i == j){return;}
  double rx,ry,rz,dd,d,ddot,mag,fx,fy,fz,nx,ny,nz,vx,vy,vz,damping,restoration,rc;
  rx = state[j*3]-state[i*3];
  ry = state[j*3+1]-state[i*3+1];
  rz = state[j*3+2]-state[i*3+2];
  dd = rx*rx+ry*ry+rz*rz;
  rc = parameters[i*2]+parameters[j*2];
  if (dd < rc*rc){
    d = std::sqrt(dd);

    nx = rx / d;
    ny = ry / d;
    nz = rz / d;

    d = rc-d;

    vx = velocities[i*3]-velocities[j*3];
    vy = velocities[i*3+1]-velocities[j*3+1];
    vz = velocities[i*3+2]-velocities[j*3+2];

    ddot = vx*nx+vy*ny+vz*nz;

    // computed per-pair from the actual colliding masses, so any
    //  number of species (Pd, Ni, H, ...) works without special-casing
    damping = ParticleSystem::damping(parameters[i*2+1],parameters[j*2+1]);
    restoration = ParticleSystem::restoration(parameters[i*2+1],parameters[j*2+1]);

    mag = -damping*ddot-restoration*d;

    // Coefficient of restitution and linear–dashpot model revisited
    //  Thomas Schwager · Thorsten Pösche
    mag = std::min(0.0,mag);

    if (mag == 0){return;}

    fx = mag*nx;
    fy = mag*ny;
    fz = mag*nz;

    forces[i*3] += fx;
    forces[i*3+1] += fy;
    forces[i*3+2] += fz;

    forces[j*3] -= fx;
    forces[j*3+1] -= fy;
    forces[j*3+2] -= fz;
  }
  else if (
    pdAttraction > 0.0 &&
    parameters[i*2+1] == mass && parameters[j*2+1] == mass
  ){
    // attractive force between big (Pd) particles, tapering linearly
    //  from pdAttraction at contact to zero at the cutoff distance
    double dOuter = rc*attractionRangeMultiplier;
    if (dd < dOuter*dOuter){
      d = std::sqrt(dd);
      nx = rx / d;
      ny = ry / d;
      nz = rz / d;

      mag = pdAttraction*(dOuter-d)/(dOuter-rc);

      fx = mag*nx;
      fy = mag*ny;
      fz = mag*nz;

      forces[i*3] += fx;
      forces[i*3+1] += fy;
      forces[i*3+2] += fz;

      forces[j*3] -= fx;
      forces[j*3+1] -= fy;
      forces[j*3+2] -= fz;
    }
  }
}

void ParticleSystem::cellCollisions(
  int64_t a1, int64_t b1, int64_t c1,
  int64_t a2, int64_t b2, int64_t c2
){
  if (
    a1 < 0 || a1 >= Nc || b1 < 0 || b1 >= Nc || c1 < 0 || c1 >= Nc ||
    a2 < 0 || a2 >= Nc || b2 < 0 || b2 >= Nc || c2 < 0 || c2 >= Nc
  ){
    return;
  }
  uint64_t cell1 = a1*Nc*Nc+b1*Nc+c1;
  uint64_t cell2 = a2*Nc*Nc+b2*Nc+c2;

  uint64_t p1 = cells[cell1];
  uint64_t p2 = cells[cell2];

  if (p1 == NULL_INDEX || p2 == NULL_INDEX){
    return;
  }

  while (p1 != NULL_INDEX){
    p2 = cells[cell2];
    while(p2 != NULL_INDEX){
        handleCollision(p1,p2);
        p2 = list[p2];
    }
    p1 = list[p1];
  }
}

double ParticleSystem::orderParameter(){
  double step = radius;
  int nl = 0;
  int ns = 0;

  for (int i = 0; i < size(); i++){
    if (parameters[i*2] == radius){
      nl += 1;
    }
    else{
      ns += 1;
    }
  }

  double mu = nl / float(nl+ns);
  double y = 0.0;
  double sigma = 0.0;
  double A = 0.0;
  while (y < Ly){
    int l = 0;
    int s = 0;
    for (int i = 0; i < size(); i++){
      if (state[i*3+1] >= y && state[i*3+1] <= y+step){
        if (parameters[i*2] == radius){
          l += 1;
        }
        else{
          s += 1;
        }
      }
    }

    if (l==0 && s==0){
      y+=step; continue;
    }

    double fi = l / float(l+s);
    double Ai = l+s;
    sigma += Ai*(fi-mu)*(fi-mu);
    A += Ai;
    y += step;
  }
  double muc = (1.0-mu)*(1.0-mu);
  return (sigma/A)/muc;
}

void ParticleSystem::applyForce(double fx, double fy, double fz){
  for (int i = 0; i < size(); i++){
    forces[i*3] += fx;
    forces[i*3+1] += fy;
    forces[i*3+2] += fz;
  }
}

void ParticleSystem::setCoeffientOfRestitution(double c){

  coefficientOfRestitution = c;
}

void ParticleSystem::newTimeStepStates(double oldDt, double newDt){
  for (int i = 0; i < size(); i++){
    for (int k = 0; k < 3; k++){
      double delta = state[i*3+k]-lastState[i*3+k];
      lastState[i*3+k] = state[i*3+k]-(newDt/oldDt)*delta;
    }
  }
}

void ParticleSystem::step(){
  clock_t tic = clock();
  resetLists();
  populateLists();
  float setup = (clock()-tic)/float(CLOCKS_PER_SEC);
  tic = clock();

  // half-neighbourhood cell stencil: self plus the 13 "forward" of the
  //  26 neighbouring cells in 3D, so every pair of nearby cells is
  //  checked exactly once as (a,b,c) sweeps the whole grid
  for (int a = 0; a < Nc; a++){
    for (int b = 0; b < Nc; b++){
      for (int c = 0; c < Nc; c++){

        cellCollisions(a,b,c, a,b,c);

        cellCollisions(a,b,c, a+1,b,c);
        cellCollisions(a,b,c, a+1,b,c+1);
        cellCollisions(a,b,c, a+1,b,c-1);
        cellCollisions(a,b,c, a+1,b+1,c);
        cellCollisions(a,b,c, a+1,b+1,c+1);
        cellCollisions(a,b,c, a+1,b+1,c-1);
        cellCollisions(a,b,c, a+1,b-1,c);
        cellCollisions(a,b,c, a+1,b-1,c+1);
        cellCollisions(a,b,c, a+1,b-1,c-1);
        cellCollisions(a,b,c, a,b+1,c);
        cellCollisions(a,b,c, a,b+1,c+1);
        cellCollisions(a,b,c, a,b+1,c-1);
        cellCollisions(a,b,c, a,b,c+1);

      }
    }
  }

  float col = (clock()-tic)/float(CLOCKS_PER_SEC);
  tic = clock();

  double cc = drag*dt/2.0;

  double damping, restoration;

  shakerDisplacement += (2.0*M_PI/shakerPeriod)*shakerAmplitude*std::sin(2.0*M_PI*shakerTime/shakerPeriod)*dt;

  for (int i = 0; i < size(); i++){

    damping = ParticleSystem::damping(parameters[i*2+1],parameters[i*2+1]);
    restoration = ParticleSystem::restoration(parameters[i*2+1],parameters[i*2+1]);

    double ct = cc/parameters[i*2+1];
    double bt = 1.0 / (1.0 + ct);
    double at = (1.0-ct)*bt;

    double x = state[i*3];
    double y = state[i*3+1];
    double z = state[i*3+2];

    double xp = lastState[i*3];
    double yp = lastState[i*3+1];
    double zp = lastState[i*3+2];

    if (y - parameters[2*i] <= (shakerDisplacement + shakerAmplitude)){
      double mag = (shakerDisplacement + shakerAmplitude)+parameters[2*i]-y;
      double f = 10*std::abs(mag)*restoration - damping*velocities[i*3+1];
      forces[i*3+1] += f;
    }

    double ax = forces[i*3];
    double ay = forces[i*3+1]-9.81*parameters[i*2+1];
    double az = forces[i*3+2];

    state[i*3] = 2.0*bt*x - at*xp + (bt*dtdt/parameters[i*2+1])*ax;
    state[i*3+1] = 2.0*bt*y - at*yp + (bt*dtdt/parameters[i*2+1])*ay;
    state[i*3+2] = 2.0*bt*z - at*zp + (bt*dtdt/parameters[i*2+1])*az;

    lastState[i*3] = x;
    lastState[i*3+1] = y;
    lastState[i*3+2] = z;

    double vx = state[i*3]-lastState[i*3];
    double vy = state[i*3+1]-lastState[i*3+1];
    double vz = state[i*3+2]-lastState[i*3+2];

    velocities[i*3] = vx/dt;
    velocities[i*3+1] = vy/dt;
    velocities[i*3+2] = vz/dt;

    double ux = 0.0; double uy = 0.0; double uz = 0.0;
    double newX = state[i*3]; double newY = state[i*3+1]; double newZ = state[i*3+2];
    bool flag = false;

    // kill the particle's motion if it's outside the box, on any wall
    if (state[i*3]-parameters[2*i] < 0 || state[i*3]+parameters[2*i] > Lx){
      ux = -0.9*vx;
      newX = (state[i*3]-parameters[2*i] < 0) ? parameters[2*i] : Lx-parameters[2*i];
      flag = true;
    }

    if (state[i*3+2]-parameters[2*i] < 0 || state[i*3+2]+parameters[2*i] > Lz){
      uz = -0.9*vz;
      newZ = (state[i*3+2]-parameters[2*i] < 0) ? parameters[2*i] : Lz-parameters[2*i];
      flag = true;
    }

    // the floor is handled by the shaker spring above; only the ceiling
    //  needs a hard reflection
    if (state[i*3+1]+parameters[2*i] > Ly){
      uy = -0.5*vy;
      newY = Ly-parameters[2*i];
      flag = true;
    }

    if (flag){
      state[i*3] = newX+ux;
      state[i*3+1] = newY+uy;
      state[i*3+2] = newZ+uz;
    }
  }

  for (int i = 0; i < size(); i++){
    forces[i*3] = 0.0;
    forces[i*3+1] = 0.0;
    forces[i*3+2] = 0.0;

    for (int k = 0; k < 3; k++){floatState[i*4+k] = float(state[i*3+k]);}
    floatState[i*4+3] = float(parameters[i*2]);
  }

  shakerTime += dt;

  double updates = (clock()-tic)/double(CLOCKS_PER_SEC);
  tic = clock();
}

void ParticleSystem::addParticle(){
  // Pd/Ni count only; grows toward the nParticles cap, inserted just
  //  before the trailing H block so H particles are left untouched
  int nPdNi = size()-nH;

  if (nPdNi == nParticles){return;}

  double x = U(generator)*(Lx-2*radius)+radius;
  double y = U(generator)*(Ly-2*radius)+radius;
  double z = U(generator)*(Lz-2*radius)+radius;

  double r = nPdNi%2 == 0 ? radius : radius / radiusRatio;
  double m = nPdNi%2 == 0 ? mass : mass / massRatio;

  insertParticle(nPdNi,x,y,z,r,m);
}

void ParticleSystem::setHCount(uint64_t target){
  while (nH < target && size() < maxCapacity){
    double x = U(generator)*(Lx-2*radius)+radius;
    double y = U(generator)*(Ly-2*radius)+radius;
    double z = U(generator)*(Lz-2*radius)+radius;
    addParticle(x,y,z,hRadius,hMass);
    nH++;
  }
  while (nH > target){
    removeParticle(size()-1);
  }
}

void ParticleSystem::randomise(double propBig){
  // only the Pd/Ni block (everything before the trailing H particles)
  //  is reassigned; H particles are untouched
  int nPdNi = size()-nH;
  int nBig = std::floor(propBig*nPdNi);
  int nSmall = nPdNi-nBig;

  for (int i = 0; i < nPdNi; i++){

    bool coin = U(generator) > 0.5;

    if (nSmall == 0 && nBig > 0){
      parameters[i*2] = radius;
      parameters[i*2+1] = mass;
      nBig--;
    }
    else if (nBig > 0 && coin){
      parameters[i*2] = radius;
      parameters[i*2+1] = mass;
      nBig--;
    }
    else if (nSmall > 0){
      parameters[i*2] = radius/radiusRatio;
      parameters[i*2+1] = mass/massRatio;
      nSmall--;
    }

    floatState[i*4+3] = parameters[i*2];

  }
}

void ParticleSystem::addParticle(
  double x,
  double y,
  double z,
  double r,
  double m
){

  int i = size();

  state.push_back(x);
  state.push_back(y);
  state.push_back(z);

  floatState[i*4] = x;
  floatState[i*4+1] = y;
  floatState[i*4+2] = z;
  floatState[i*4+3] = r;

  lastState.push_back(x);
  lastState.push_back(y);
  lastState.push_back(z);

  parameters.push_back(r);
  parameters.push_back(m);

  forces.push_back(0.0);
  forces.push_back(0.0);
  forces.push_back(0.0);

  velocities.push_back(0.0);
  velocities.push_back(0.0);
  velocities.push_back(0.0);
}

void ParticleSystem::insertParticle(
  uint64_t idx,
  double x,
  double y,
  double z,
  double r,
  double m
){

  state.insert(state.begin()+3*idx, {x,y,z});
  lastState.insert(lastState.begin()+3*idx, {x,y,z});
  parameters.insert(parameters.begin()+2*idx, {r,m});
  forces.insert(forces.begin()+3*idx, {0.0,0.0,0.0});
  velocities.insert(velocities.begin()+3*idx, {0.0,0.0,0.0});

  // shift the trailing floatState entries (the fixed C-array has no
  //  insert of its own) up by one slot to make room at idx
  int n = size();
  for (int i = n-1; i > idx; i--){
    floatState[i*4]   = floatState[(i-1)*4];
    floatState[i*4+1] = floatState[(i-1)*4+1];
    floatState[i*4+2] = floatState[(i-1)*4+2];
    floatState[i*4+3] = floatState[(i-1)*4+3];
  }
  floatState[idx*4] = x;
  floatState[idx*4+1] = y;
  floatState[idx*4+2] = z;
  floatState[idx*4+3] = r;
}

void ParticleSystem::removeParticle(uint64_t i){
  if (state.size() >= 3*i){
    // an index in the trailing block beyond nParticles is a proteium
    //  (H) particle; keep the H count in sync however removal happens
    if (i >= size()-nH){nH--;}

    state.erase(
      state.begin()+3*i,
      state.begin()+3*i+3
    );

    lastState.erase(
      lastState.begin()+3*i,
      lastState.begin()+3*i+3
    );

    parameters.erase(
      parameters.begin()+2*i,
      parameters.begin()+2*i+2
    );

    forces.erase(
      forces.begin()+3*i,
      forces.begin()+3*i+3
    );

    velocities.erase(
      velocities.begin()+3*i,
      velocities.begin()+3*i+3
    );
  }
}

void ParticleSystem::oneBigOnBottom(){
  parameters[0] = radius;
  parameters[1] = mass;
  floatState[0*4+3] = radius;
  // leave the trailing proteium (H) block untouched
  int nPdNi = size()-nH;
  for (int i = 1; i < nPdNi; i++){
    parameters[i*2] = radius/radiusRatio;
    parameters[i*2+1] = mass/massRatio;
    floatState[i*4+3] = parameters[i*2];
  }

  state[0] = Lx/2.0;
  state[1] = 2*radius;
  state[2] = Lz/2.0;

  double r = 2.0*radius/radiusRatio;
  int n = std::floor(Lx/r);
  int j = 0;
  double y = r+3*radius;
  for (int i = 1; i < nPdNi; i++){
    if (j < n){
      state[i*3] = j*r+radius/2.;
      state[i*3+1] = y;
      state[i*3+2] = Lz/2.0;
      j++;
    }
    else{
      j = 0;
      y += r*1.1;
      state[i*3] = j*r+radius/2.;
      state[i*3+1] = y;
      state[i*3+2] = Lz/2.0;
    }
  }

  for (int i = 0; i < nPdNi; i++){
    lastState[i*3] = state[i*3];
    lastState[i*3+1] = state[i*3+1];
    lastState[i*3+2] = state[i*3+2];
    floatState[i*4] = state[i*3];
    floatState[i*4+1] = state[i*3+1];
    floatState[i*4+2] = state[i*3+2];
  }

  shakerDisplacement = 0;
  shakerPeriod = 1.0;
  shakerAmplitude = 0;
}

void ParticleSystem::one(){

  while (size() > 0){
    removeParticle();
  }

  addParticle();

  parameters[0] = radius;
  parameters[1] = mass;
  floatState[0*4+3] = radius;

  state[0] = Lx/2.0;
  state[1] = Ly-2*radius;
  state[2] = Lz/2.0;

  lastState[0] = state[0];
  lastState[1] = state[1];
  lastState[2] = state[2];
  floatState[0] = state[0];
  floatState[1] = state[1];
  floatState[2] = state[2];

  shakerDisplacement = 0;
  shakerPeriod = 1.0;
  shakerAmplitude = 0;

}
