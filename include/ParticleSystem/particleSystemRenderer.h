#ifndef PARTICLESYSTEMRENDERER_H
#define PARTICLESYSTEMRENDERER_H

#include <glm/glm.hpp>

class ParticleSystemRenderer {
public:
  ParticleSystemRenderer(int sizeHint, int trackLength = 60*100)
  : nParticles(sizeHint)
  {

    trackedParticle = NULL_INDEX;
    // x,y,z,fade-time per trail point
    track = new float [trackLength*4];
    this->trackLength = trackLength;

    initialiseGL();
  }

  void setProjection(glm::mat4 p);
  void draw(ParticleSystem & p, uint64_t frameId, float resX, float resY);

  ~ParticleSystemRenderer(){
    // kill some GL stuff
    glDeleteProgram(particleShader);
    glDeleteProgram(shakerShader);
    glDeleteProgram(trackShader);

    glDeleteBuffers(1,&offsetVBO);
    glDeleteBuffers(1,&sphereVertVBO);
    glDeleteBuffers(1,&sphereNormVBO);
    glDeleteBuffers(1,&sphereEBO);
    glDeleteBuffers(1,&shakerVBO);
    glDeleteBuffers(1,&trackVBO);

    glDeleteVertexArrays(1,&vertVAO);
    glDeleteVertexArrays(1,&shakerVAO);
    glDeleteVertexArrays(1,&trackVAO);

    free(track);
  }

  // rayOrigin/rayDir define a world-space pick ray (see Camera3D::screenRay)
  void click(ParticleSystem & p, glm::vec3 rayOrigin, glm::vec3 rayDir);
  void beginTracking(ParticleSystem & p, uint64_t i);
  void updatedTrack(ParticleSystem & p);

private:
  int nParticles;
  // GL data members
  float * floatState;
  GLuint particleShader, offsetVBO, vertVAO;
  glm::mat4 projection;

  // instanced low-poly sphere mesh, generated once in initialiseGL()
  GLuint sphereVertVBO, sphereNormVBO, sphereEBO;
  int sphereIndexCount;

  GLuint shakerShader, shakerVAO, shakerVBO;
  // horizontal quad (x,z pairs), rewritten each frame to match the box size
  float shakerVertices[6*2];

  float * track;
  int trackLength;
  uint64_t trackedParticle;

  GLuint trackShader, trackVAO, trackVBO;

  void initialiseGL();
  void buildSphereMesh(int stacks, int slices);

};

#endif
