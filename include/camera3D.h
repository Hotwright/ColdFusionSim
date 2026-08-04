#ifndef CAMERA3D_H
#define CAMERA3D_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

/*
  A perspective orbit camera: orbits a target point at a given
  distance, controlled by yaw/pitch angles.
*/
class Camera3D {
public:
  Camera3D(
    int resX, int resY,
    glm::vec3 target = glm::vec3(0.25f,0.5f,0.25f)
  )
  : resolution(resX,resY), target(target), distance(1.6f),
    yaw(45.0f), pitch(20.0f), fov(45.0f),
    nearPlane(0.01f), farPlane(100.0f)
  {
    update();
  }

  glm::mat4 getVP(){return vp;}
  glm::mat4 getView(){return view;}
  glm::mat4 getProjection(){return projection;}
  glm::vec3 getEye(){return eye;}
  glm::vec3 getTarget(){return target;}
  float getYaw(){return yaw;}
  float getPitch(){return pitch;}
  float getDistance(){return distance;}

  void orbit(float dYaw, float dPitch){
    yaw += dYaw;
    pitch += dPitch;
    if (pitch > 89.0f){pitch = 89.0f;}
    if (pitch < -89.0f){pitch = -89.0f;}
    update();
  }

  void zoom(float dz){
    distance -= dz*0.15f*distance;
    if (distance < 0.05f){distance = 0.05f;}
    if (distance > 20.0f){distance = 20.0f;}
    update();
  }

  void setTarget(glm::vec3 t){target = t; update();}

  // slides the orbit target across the view plane (right/up in camera
  //  space), so zooming in afterward focuses on wherever was panned to;
  //  scaled by distance so a drag feels the same size on screen at any zoom
  void pan(float dx, float dy){
    glm::vec3 right(view[0][0],view[1][0],view[2][0]);
    glm::vec3 up(view[0][1],view[1][1],view[2][1]);
    float scale = distance*0.0015f;
    target += -right*dx*scale+up*dy*scale;
    update();
  }

  void reset(glm::vec3 t, float d, float y, float p){
    target = t; distance = d; yaw = y; pitch = p;
    update();
  }

  // unproject a screen pixel (SFML coords, y down) into a world-space ray
  void screenRay(float x, float y, glm::vec3 & outOrigin, glm::vec3 & outDir){
    float ndcX = 2.0f*x/resolution.x-1.0f;
    float ndcY = 1.0f-2.0f*y/resolution.y;

    glm::mat4 invVP = glm::inverse(vp);

    glm::vec4 nearP = invVP*glm::vec4(ndcX,ndcY,-1.0f,1.0f);
    glm::vec4 farP  = invVP*glm::vec4(ndcX,ndcY, 1.0f,1.0f);
    nearP /= nearP.w;
    farP /= farP.w;

    outOrigin = glm::vec3(nearP);
    outDir = glm::normalize(glm::vec3(farP-nearP));
  }

private:

  void update(){
    float ry = glm::radians(yaw);
    float rp = glm::radians(pitch);

    glm::vec3 offset(
      distance*cos(rp)*cos(ry),
      distance*sin(rp),
      distance*cos(rp)*sin(ry)
    );

    eye = target+offset;

    view = glm::lookAt(eye,target,glm::vec3(0.0f,1.0f,0.0f));
    projection = glm::perspective(
      glm::radians(fov),
      resolution.x/resolution.y,
      nearPlane,
      farPlane
    );

    vp = projection*view;
  }

  glm::vec2 resolution;
  glm::vec3 target;
  glm::vec3 eye;

  float distance;
  float yaw;
  float pitch;
  float fov;
  float nearPlane;
  float farPlane;

  glm::mat4 view;
  glm::mat4 projection;
  glm::mat4 vp;
};

#endif
