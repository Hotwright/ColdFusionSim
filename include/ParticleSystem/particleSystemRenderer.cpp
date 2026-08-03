#include <ParticleSystem/particleSystemRenderer.h>
#include <vector>
#include <cmath>

void ParticleSystemRenderer::setProjection(glm::mat4 p){
  projection = p;
  glUseProgram(particleShader);
  glUniformMatrix4fv(
    glGetUniformLocation(particleShader,"proj"),
    1,
    GL_FALSE,
    &projection[0][0]
  );

  glUseProgram(shakerShader);
  glUniformMatrix4fv(
    glGetUniformLocation(shakerShader,"proj"),
    1,
    GL_FALSE,
    &projection[0][0]
  );

  glUseProgram(trackShader);

  glUniformMatrix4fv(
    glGetUniformLocation(trackShader,"proj"),
    1,
    GL_FALSE,
    &projection[0][0]
  );

}

void ParticleSystemRenderer::buildSphereMesh(int stacks, int slices){
  std::vector<float> verts;
  std::vector<float> norms;
  std::vector<unsigned int> indices;

  for (int i = 0; i <= stacks; i++){
    double phi = M_PI*double(i)/double(stacks);
    double y = std::cos(phi);
    double r = std::sin(phi);
    for (int j = 0; j <= slices; j++){
      double theta = 2.0*M_PI*double(j)/double(slices);
      double x = r*std::cos(theta);
      double z = r*std::sin(theta);
      verts.push_back(float(x));
      verts.push_back(float(y));
      verts.push_back(float(z));
      norms.push_back(float(x));
      norms.push_back(float(y));
      norms.push_back(float(z));
    }
  }

  for (int i = 0; i < stacks; i++){
    for (int j = 0; j < slices; j++){
      unsigned int a = i*(slices+1)+j;
      unsigned int b = a+slices+1;
      indices.push_back(a);
      indices.push_back(b);
      indices.push_back(a+1);
      indices.push_back(a+1);
      indices.push_back(b);
      indices.push_back(b+1);
    }
  }

  sphereIndexCount = indices.size();

  glGenBuffers(1,&sphereVertVBO);
  glBindBuffer(GL_ARRAY_BUFFER,sphereVertVBO);
  glBufferData(GL_ARRAY_BUFFER,sizeof(float)*verts.size(),verts.data(),GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER,0);

  glGenBuffers(1,&sphereNormVBO);
  glBindBuffer(GL_ARRAY_BUFFER,sphereNormVBO);
  glBufferData(GL_ARRAY_BUFFER,sizeof(float)*norms.size(),norms.data(),GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER,0);

  glGenBuffers(1,&sphereEBO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,sphereEBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(unsigned int)*indices.size(),indices.data(),GL_STATIC_DRAW);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
}

void ParticleSystemRenderer::initialiseGL(){
  // a buffer of particle states: x,y,z,radius per instance
  glGenBuffers(1,&offsetVBO);
  glBindBuffer(GL_ARRAY_BUFFER,offsetVBO);
  glBufferData(GL_ARRAY_BUFFER,sizeof(float)*nParticles*4,NULL,GL_DYNAMIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER,0);

  buildSphereMesh(8,12);

  // setup an array object: base sphere mesh + per-instance offset
  glGenVertexArrays(1,&vertVAO);
  glBindVertexArray(vertVAO);

  glBindBuffer(GL_ARRAY_BUFFER,sphereVertVBO);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,3*sizeof(float),(void*)0);

  glEnableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER, offsetVBO);
  glVertexAttribPointer(1,4,GL_FLOAT,GL_FALSE,4*sizeof(float),(void*)0);
  glVertexAttribDivisor(1,1);

  glBindBuffer(GL_ARRAY_BUFFER,sphereNormVBO);
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2,3,GL_FLOAT,GL_FALSE,3*sizeof(float),(void*)0);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,sphereEBO);

  glBindBuffer(GL_ARRAY_BUFFER,0);
  glBindVertexArray(0);

  glError("initialised particles");

  particleShader = glCreateProgram();
  compileShader(particleShader,particleVertexShader,particleFragmentShader);
  glUseProgram(particleShader);

  glUniformMatrix4fv(
    glGetUniformLocation(particleShader,"proj"),
    1,
    GL_FALSE,
    &projection[0][0]
  );

  // shaker: a horizontal quad (x,z pairs), sized to the box in draw()

  shakerShader = glCreateProgram();
  compileShader(shakerShader,shakerVertexShader,shakerFragmentShader);
  glUseProgram(shakerShader);

  glUniformMatrix4fv(
    glGetUniformLocation(shakerShader,"proj"),
    1,
    GL_FALSE,
    &projection[0][0]
  );

  glUniform4f(
    glGetUniformLocation(shakerShader,"u_colour"),
    0.35f,0.35f,0.35f,0.55f
  );

  for (int i = 0; i < 12; i++){shakerVertices[i] = 0.0f;}

  glGenVertexArrays(1,&shakerVAO);
  glGenBuffers(1,&shakerVBO);
  glBindVertexArray(shakerVAO);
  glBindBuffer(GL_ARRAY_BUFFER,shakerVBO);
  glBufferData(GL_ARRAY_BUFFER,sizeof(float)*6*2,shakerVertices,GL_DYNAMIC_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,2*sizeof(float),0);
  glBindBuffer(GL_ARRAY_BUFFER,0);
  glBindVertexArray(0);

  glError();
  glBufferStatus();

  // track: x,y,z,fade-time per point

  glGenVertexArrays(1,&trackVAO);
  glGenBuffers(1,&trackVBO);
  glBindVertexArray(trackVAO);
  glBindBuffer(GL_ARRAY_BUFFER,trackVBO);
  glBufferData(GL_ARRAY_BUFFER,sizeof(float)*trackLength*4,track,GL_DYNAMIC_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,4*sizeof(float),0);
  glBindBuffer(GL_ARRAY_BUFFER,0);
  glBindVertexArray(0);

  trackShader = glCreateProgram();
  compileShader(trackShader,trackVertexShader,trackFragmentShader);
  glUseProgram(trackShader);

  glUniformMatrix4fv(
    glGetUniformLocation(trackShader,"proj"),
    1,
    GL_FALSE,
    &projection[0][0]
  );

  glUniform1f(
    glGetUniformLocation(trackShader,"time"),
    trackLength
  );

}

void ParticleSystemRenderer::draw(
  ParticleSystem & p,
  uint64_t frameId,
  float resX,
  float resY
){
  glUseProgram(particleShader);

  glUniform1i(
    glGetUniformLocation(particleShader,"tracked"),
    trackedParticle != NULL_INDEX ? trackedParticle : -1
  );

  // midpoints between adjacent species' radii, used to colour particles:
  //  red = big/Pd, blue = small/Ni, green = tiny/H (proteium)
  glUniform1f(
    glGetUniformLocation(particleShader,"sizeThresholdHi"),
    0.5*(p.radius+p.radius/p.radiusRatio)
  );

  glUniform1f(
    glGetUniformLocation(particleShader,"sizeThresholdLo"),
    0.5*(p.radius/p.radiusRatio+p.hRadius)
  );

  glBindBuffer(GL_ARRAY_BUFFER,offsetVBO);
  glBufferSubData(GL_ARRAY_BUFFER,0,sizeof(float)*nParticles*4,&p.floatState[0]);
  glBindBuffer(GL_ARRAY_BUFFER,0);

  glError("particles buffers");

  glBindVertexArray(vertVAO);
  glDrawElementsInstanced(GL_TRIANGLES,sphereIndexCount,GL_UNSIGNED_INT,0,p.size());
  glBindVertexArray(0);

  glError("draw particles");

  glUseProgram(shakerShader);

  glUniform1f(
    glGetUniformLocation(shakerShader,"offset"),
    p.shakerDisplacement+p.shakerAmplitude
  );

  // rewrite the floor quad to match the box footprint each frame
  //  (cheap: 6 vertices) since Lx/Lz live on the particle system
  shakerVertices[0] = 0.0f;      shakerVertices[1] = 0.0f;
  shakerVertices[2] = p.Lx;      shakerVertices[3] = 0.0f;
  shakerVertices[4] = p.Lx;      shakerVertices[5] = p.Lz;
  shakerVertices[6] = 0.0f;      shakerVertices[7] = 0.0f;
  shakerVertices[8] = p.Lx;      shakerVertices[9] = p.Lz;
  shakerVertices[10] = 0.0f;     shakerVertices[11] = p.Lz;

  glBindBuffer(GL_ARRAY_BUFFER,shakerVBO);
  glBufferSubData(GL_ARRAY_BUFFER,0,sizeof(float)*6*2,shakerVertices);
  glBindBuffer(GL_ARRAY_BUFFER,0);

  glBindVertexArray(shakerVAO);
  glDrawArrays(GL_TRIANGLES,0,6);
  glBindVertexArray(0);

  glError("draw shaker");

  if (trackedParticle != NULL_INDEX){
    glUseProgram(trackShader);

    updatedTrack(p);

    glBindBuffer(GL_ARRAY_BUFFER,trackVBO);
    glBufferSubData(GL_ARRAY_BUFFER,0,sizeof(float)*trackLength*4,&track[0]);
    glBindBuffer(GL_ARRAY_BUFFER,0);

    glBindVertexArray(trackVAO);
    glDrawArrays(GL_POINTS,0,trackLength);
    glBindVertexArray(0);

  }

}

void ParticleSystemRenderer::updatedTrack(ParticleSystem & p){
  if (trackedParticle != NULL_INDEX){

    float oldx = track[0]; float oldy = track[1]; float oldz = track[2];
    for (int t = 1; t < trackLength-1; t++){
      float x = track[t*4]; float y = track[t*4+1]; float z = track[t*4+2];
      track[t*4] = oldx;
      track[t*4+1] = oldy;
      track[t*4+2] = oldz;
      oldx = x; oldy = y; oldz = z;
    }

    track[0] = p.floatState[trackedParticle*4];
    track[1] = p.floatState[trackedParticle*4+1];
    track[2] = p.floatState[trackedParticle*4+2];

  }
}

void ParticleSystemRenderer::beginTracking(ParticleSystem & p, uint64_t i){
  if (i == trackedParticle){return;}

  trackedParticle = i;
  for (int t = 0; t < trackLength; t++){
    track[t*4] = p.floatState[trackedParticle*4];
    track[t*4+1] = p.floatState[trackedParticle*4+1];
    track[t*4+2] = p.floatState[trackedParticle*4+2];
    track[t*4+3] = trackLength-t;
  }

}

void ParticleSystemRenderer::click(ParticleSystem & p, glm::vec3 rayOrigin, glm::vec3 rayDir){

  uint64_t best = NULL_INDEX;
  double bestT = 1e18;

  for (int i = 0; i < p.size(); i++){
    glm::vec3 c(p.state[i*3],p.state[i*3+1],p.state[i*3+2]);
    double r = p.parameters[i*2];

    glm::vec3 oc = rayOrigin-c;
    double b = glm::dot(oc,rayDir);
    double cterm = glm::dot(oc,oc)-r*r;
    double disc = b*b-cterm;

    if (disc < 0.0){continue;}

    double t = -b-std::sqrt(disc);
    if (t > 0.0 && t < bestT){
      bestT = t;
      best = i;
    }
  }

  if (best != NULL_INDEX){
    return beginTracking(p,best);
  }

  trackedParticle = NULL_INDEX;
}
