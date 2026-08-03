#ifndef SHADERS_H
#define SHADERS_H

const char * trackVertexShader = "#version 330 core\n"
  "layout(location=0) in vec4 a_position;\n"
  "uniform mat4 proj; uniform float time;\n"
  "out float alpha;\n"
  "void main(){\n"
    "gl_Position = proj*vec4(a_position.xyz,1.0);\n"
    "gl_PointSize = 4.0;\n"
    "alpha = a_position.w/time;"
  "}";

const char * trackFragmentShader = "#version 330 core\n"
  "out vec4 colour; in float alpha;\n"
  "void main(){"
  " vec2 c = 2.0*gl_PointCoord-1.0;\n"
  " float d = length(c);\n"
  // bit of simple AA
  " float a = 1.0-smoothstep(0.99,1.01,d);\n"
  " colour = vec4(0.0,0.0,0.0,alpha);\n"
  " if (colour.a == 0.0){discard;}"
  "}";

// a horizontal quad (x,z) at y=offset, spanning the box footprint
const char * shakerVertexShader = "#version 330 core\n"
  "layout(location=0) in vec2 a_position;\n"
  "\n"
  "uniform mat4 proj;\n"
  "uniform float offset;\n"
  "void main(){\n"
  " gl_Position = proj*vec4(a_position.x,offset,a_position.y,1.0);"
  " \n"
  "}";

const char * shakerFragmentShader = "#version 330 core\n"
  "uniform vec4 u_colour; out vec4 colour;\n"
  "void main(){colour=u_colour;\n}";

const char * boxVertexShader = "#version 330 core\n"
  "layout(location=0) in vec2 a_position;\n"
  "uniform vec3 colour; out vec4 o_colour;\n"
  "uniform mat4 proj;\n"
  "void main(){\n"
  " vec4 pos = vec4(a_position.xy,0.0,1.0);\n"
  " o_colour = vec4(colour,.5);\n"
  " gl_Position = pos;\n"
  "}";

const char * boxFragmentShader = "#version 330 core\n"
  "in vec4 o_colour; out vec4 colour;\n"
  "void main(){colour=o_colour;\n}";

const char * sliderVertexShader = "#version 330 core\n"
  "layout(location=0) in vec2 a_position;\n"
  "uniform mat4 proj;\n"
  "void main(){\n"
  " gl_Position = proj*vec4(a_position.xy,0.0,1.0);\n"
  "}";

const char * sliderFragmentShader = "#version 330 core\n"
  "uniform vec4 frameColour; uniform vec4 fillColour;\n"
  "uniform vec2 state;\n"
  "out vec4 colour;\n"
  "void main(){"
    "if (gl_FragCoord.x > state.x && gl_FragCoord.x < state.y){colour=fillColour;}"
    "else {colour = frameColour;}"
  "}";

const char * buttonVertexShader = "#version 330 core\n"
  "layout(location=0) in vec2 a_position;\n"
  "uniform mat4 proj;\n"
  "void main(){\n"
  " gl_Position = proj*vec4(a_position.xy,0.0,1.0);\n"
  "}";

const char * buttonFragmentShader = "#version 330 core\n"
  "uniform vec4 frameColour; uniform vec4 fillColour;\n"
  "uniform int state; uniform float alpha;\n"
  "out vec4 colour;\n"
  "void main(){"
    "colour=fillColour*alpha+frameColour;\n"
  "}";


// instanced, lit low-poly sphere impostor: a_position/a_normal come from
//  one shared unit-sphere mesh, a_offset carries each particle's
//  (x,y,z,radius) so the same mesh is reused for every instance
const char * particleVertexShader = "#version 330 core\n"
  "precision highp float;\n"
  "layout(location = 0) in vec3 a_position;\n"
  "layout(location = 1) in vec4 a_offset;\n"
  "layout(location = 2) in vec3 a_normal;\n"
  "uniform mat4 proj;\n"
  "uniform int tracked;\n"
  "uniform float sizeThresholdHi;\n"
  "uniform float sizeThresholdLo;\n"
  "out vec3 o_normal;\n"
  "out vec4 o_colour;\n"
  "void main(){\n"
  " vec3 worldPos = a_position*a_offset.w+a_offset.xyz;\n"
  " gl_Position = proj*vec4(worldPos,1.0);\n"
  " o_normal = a_normal;\n"
  "if (tracked == gl_InstanceID){o_colour = vec4(1.0,1.0,0.0,1.0);}"
  "else if (a_offset.w > sizeThresholdHi){ o_colour = vec4(1.0,0.0,0.0,1.0); }" // Pd, red
  "else if (a_offset.w > sizeThresholdLo){ o_colour = vec4(0.0,0.0,1.0,1.0); }" // Ni, blue
  "else{ o_colour = vec4(0.0,0.7,0.0,1.0); }\n" // H, green
  "}";
const char * particleFragmentShader = "#version 330 core\n"
  "in vec3 o_normal; in vec4 o_colour; out vec4 colour;\n"
  "void main(){\n"
  " vec3 n = normalize(o_normal);\n"
  " vec3 lightDir = normalize(vec3(0.4,0.8,0.5));\n"
  " float diff = max(dot(n,lightDir),0.0);\n"
  " float ambient = 0.35;\n"
  " vec3 rgb = o_colour.rgb*(ambient+diff*0.65);\n"
  " colour = vec4(rgb,o_colour.a);\n"
  "}";

#endif
