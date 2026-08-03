#include <iostream>
#include <sstream>
#include <iomanip>

#if WINDOWS
  #include <glew.c>
#else
  #include <glew.h>
#endif

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>

#include <orthoCam.h>
#include <camera3D.h>

#include <glUtils.h>
#include <utils.h>
#include <shaders.h>

#include <ParticleSystem/particleSystem.cpp>
#include <ParticleSystem/particleSystemRenderer.cpp>
#include <ParticleSystem/trajectory.cpp>
#include <Text/textRenderer.cpp>
#include <Text/popup.cpp>

#include <Menu/button.cpp>
#include <Menu/slider.cpp>

#include <time.h>
#include <random>
#include <iostream>
#include <math.h>
#include <vector>

const int resX = 1000;
const int resY = 1000;
// the 3D scene renders only into this region of the window: left of the
//  control column, and dropped down from the top by simTopMargin
const int panelWidth = 340;
const int simWidth = resX-panelWidth;
const int simTopMargin = 320;
const int simHeight = resY-simTopMargin;

// empirical atomic radii (picometres), Clementi et al. 1967
const double HYDROGEN_ATOMIC_RADIUS = 25.0;
const double NICKEL_ATOMIC_RADIUS = 124.0;
const double PALLADIUM_ATOMIC_RADIUS = 137.0;
const double PD_NI_RADIUS_RATIO = PALLADIUM_ATOMIC_RADIUS / NICKEL_ATOMIC_RADIUS;
const double PD_H_RADIUS_RATIO = PALLADIUM_ATOMIC_RADIUS / HYDROGEN_ATOMIC_RADIUS;

// standard atomic weights (u), IUPAC
const double NICKEL_ATOMIC_MASS = 58.6934;
const double PALLADIUM_ATOMIC_MASS = 106.42;
const double PD_NI_MASS_RATIO = PALLADIUM_ATOMIC_MASS / NICKEL_ATOMIC_MASS;

// the simulation box (world units)
const double BOX_LX = 0.5;
const double BOX_LY = 1.0;
const double BOX_LZ = 0.5;

const int subSamples = 60;
const float dt = (1.0 / 60.0) / subSamples;

const int saveFrequency = 1;

const int N = 1024;
// proteium (H) particles are added on top of N (the Pd/Ni cap), up to
//  this many per Pd particle, so storage must be reserved for the worst case
const double MAX_H_RATIO = 10.0;
const uint64_t CAPACITY = N + uint64_t(N*MAX_H_RATIO);
// motion parameters

// for smoothing delta numbers
uint8_t frameId = 0;
double deltas[60];
double physDeltas[60];
double renderDeltas[60];

float speed = 1.0;

void speedPopUp(Popup & popups){
  popups.clear("speed");
  std::string pos = fixedLengthNumber(std::ceil(100.0*speed),3);
  popups.post(
    FadingText(
      "Speed: "+pos+" %",
      3.0,
      resX-256.,
      resY-64*9,
      glm::vec3(0.0,0.0,0.0),
      "speed"
    )
  );
}

int main(){

  for (int i = 0; i < 60; i++){deltas[i] = 0.0;}

  sf::ContextSettings contextSettings;
  contextSettings.depthBits = 24;
  contextSettings.antialiasingLevel = 0;

  sf::RenderWindow window(
    sf::VideoMode(resX,resY),
    "Cold Fusion Sim",
    sf::Style::Close|sf::Style::Titlebar,
    contextSettings
  );
  window.setVerticalSyncEnabled(true);
  window.setFramerateLimit(60);
  window.setActive();

  glewInit();

  uint8_t debug = 0;

  // the core simulation
  ParticleSystem particles(N,CAPACITY,dt,0.25,BOX_LX,BOX_LY,BOX_LZ);
  particles.setHRadiusRatio(PD_H_RADIUS_RATIO);
  // handles rendering - separation of concerns
  ParticleSystemRenderer pRender(CAPACITY);

  sf::Clock clock;
  sf::Clock physClock, renderClock;

  glm::mat4 textProj = glm::ortho(0.0,double(resX),0.0,double(resY));

  glEnable( GL_BLEND );
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDepthFunc(GL_LEQUAL);

  // for freetype rendering
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

  // must be initialised before so the shader is in use..?
  TextRenderer textRenderer(textProj);

  Type OD("resources/fonts/","OpenDyslexic-Regular.otf",48);

  Popup popups;

  // aspect ratio matches the viewport the 3D scene actually renders
  //  into (see simWidth/simHeight), not the full window
  Camera3D camera(simWidth,simHeight,glm::vec3(BOX_LX/2.0,BOX_LY/2.0,BOX_LZ/2.0));

  glViewport(0,0,resX,resY);

  // sliders are not beautifully handeled, should really have
  //  a widget hierachy system, but won't bother for this
  //  app

  Slider amplitudeSlider(resX-300.0,resY-64.0,128.0,16.0,"Shaker Amplitude");
  amplitudeSlider.setPosition(0.5);
  amplitudeSlider.setProjection(textProj);

  Slider shakerSlider(resX-300.0,resY-64.0*2,128.0,16.0,"Shaker Period");
  shakerSlider.setPosition(0.5);
  shakerSlider.setProjection(textProj);

  Slider particlesSlider(resX-300.0,resY-64.0*3,128.0,16.0,"Particles");
  particlesSlider.setPosition(0.5);
  particlesSlider.setProjection(textProj);

  Slider proportionBigSlider(resX-300.0,resY-64.0*4,128.0,16.0,"Pd:Ni Ratio");
  proportionBigSlider.setPosition(0.5);
  proportionBigSlider.setProjection(textProj);

  Slider restitutionSlider(resX-300.0,resY-64.0*5,128.0,16.0,"Coef. Restitution");
  restitutionSlider.setPosition(0.5);
  restitutionSlider.setProjection(textProj);

  Slider massRatioSlider(resX-300.0,resY-64.0*6,128.0,16.0,"Mass Ratio");
  // default to the real Pd:Ni atomic mass ratio (slider maps [0,1] -> [1, 1+maxMassRatio])
  massRatioSlider.setPosition((PD_NI_MASS_RATIO-1.0)/2.0);
  massRatioSlider.setProjection(textProj);

  Slider radiusRatioSlider(resX-300.0,resY-64.0*7,128.0,16.0,"Size Ratio");
  // default to the real Pd:Ni atomic radius ratio (slider maps [0,1] -> [1, 1+maxRadiusRatio])
  radiusRatioSlider.setPosition((PD_NI_RADIUS_RATIO-1.0)/4.0);
  radiusRatioSlider.setProjection(textProj);

  Slider pdAttractionSlider(resX-300.0,resY-64.0*8,128.0,16.0,"Pd Attraction");
  pdAttractionSlider.setPosition(0.0);
  pdAttractionSlider.setProjection(textProj);

  Slider proteiumSlider(resX-300.0,resY-64.0*9,128.0,16.0,"Proteium Count");
  proteiumSlider.setPosition(0.0);
  proteiumSlider.setProjection(textProj);

  // buttons live directly under all the sliders above
  Button newRecording(resX-300.0,resY-64.0*10,16.0,16.0,"Record",30);
  newRecording.setState(false);
  newRecording.setProjection(textProj);

  Button randomizeButton(resX-300.0,resY-64.0*11,16.0,16.0,"Randomize",30);
  randomizeButton.setState(false);
  randomizeButton.setProjection(textProj);

  std::vector<Slider*> sliders = {
    &amplitudeSlider, &shakerSlider, &particlesSlider, &proportionBigSlider,
    &restitutionSlider, &massRatioSlider, &radiusRatioSlider,
    &pdAttractionSlider, &proteiumSlider
  };
  // whichever slider was last clicked/dragged; arrow keys nudge this one
  Slider * activeSlider = nullptr;

  double oldMouseX = 0.0;
  double oldMouseY = 0.0;

  double mouseX = resX/2.0;
  double mouseY = resY/2.0;

  double oldOrbitX = 0.0;
  double oldOrbitY = 0.0;
  bool orbiting = false;

  bool moving = false;

  bool pause = false;

  double shakerMaxPeriod = 1.0;
  double propBig = 0.5;
  double maxAmplitude = 10.0; // measured in particle radius units!
  double maxMassRatio = 2.0;
  double maxRadiusRatio = 4.0;
  double maxPdAttraction = 10.0;

  bool isRecording = false;
  Trajectory record;

  while (window.isOpen()){

    sf::Event event;
    while (window.pollEvent(event)){

      if (event.type == sf::Event::Closed){
        if ( isRecording ){record.save();}
        window.close();
      }

      if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::O){
        particlesSlider.setPosition(1.0/float(N));
        particles.one();
        pause = true;
      }

      if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::L){
        shakerSlider.smoothChangeTo(0.133,60*5,2.0);
      }


      if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::D){
        amplitudeSlider.setSmoothChange(true,60*5,2.0);
        shakerSlider.setSmoothChange(true,60*5,2.0);
        restitutionSlider.setSmoothChange(true,60*5,2.0);
      }

      if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape){
        if ( isRecording ){record.save();}
        window.close();
      }

      if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::F1){
        debug = !debug;
      }

      if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Space){
        pause = !pause;
      }


      if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::W){
        speed *= 2.0;
        if (speed > 1){
          speed = 1;
        }
        speedPopUp(popups);
      }

      if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::S){
        speed /= 2.0;
        if (speed < 0.01){
          speed = 0.01;
        }
        speedPopUp(popups);
      }


      if (event.type == sf::Event::MouseWheelScrolled){
        mouseX = event.mouseWheelScroll.x;
        mouseY = event.mouseWheelScroll.y;
        double z = event.mouseWheelScroll.delta;

        camera.zoom(z);
      }

      // middle click resets the view
      if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Middle){
        camera.reset(glm::vec3(BOX_LX/2.0,BOX_LY/2.0,BOX_LZ/2.0),1.6f,45.0f,20.0f);
      }

      // right-drag orbits the camera around the box
      if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Right){
        oldOrbitX = event.mouseButton.x;
        oldOrbitY = event.mouseButton.y;
        orbiting = true;
      }

      if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Right){
        orbiting = false;
      }

      if (event.type == sf::Event::MouseMoved && orbiting){
        sf::Vector2i pos = sf::Mouse::getPosition(window);
        camera.orbit(float(pos.x-oldOrbitX)*0.3f, float(oldOrbitY-pos.y)*0.3f);
        oldOrbitX = pos.x;
        oldOrbitY = pos.y;
      }

      if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left){
        sf::Vector2i pos = sf::Mouse::getPosition(window);

        // whichever slider gets hit becomes the arrow-key target
        for (Slider * s : sliders){
          if (s->clicked(pos.x,resY-pos.y)){activeSlider = s;}
        }

        // buttons
        newRecording.clicked(pos.x,resY-pos.y);
        randomizeButton.clicked(pos.x,resY-pos.y);

        // cast a ray through the clicked pixel and find the nearest particle it
        //  hits; only meaningful within the 3D viewport, not the control column
        //  or the top margin, and rebased to the viewport's own coordinates
        if (pos.x < simWidth && pos.y >= simTopMargin){
          glm::vec3 rayOrigin, rayDir;
          camera.screenRay(pos.x,pos.y-simTopMargin,rayOrigin,rayDir);

          pRender.click(particles,rayOrigin,rayDir);
        }

        oldMouseX = pos.x;
        oldMouseY = pos.y;
      }

      if (event.type == sf::Event::MouseMoved && sf::Mouse::isButtonPressed(sf::Mouse::Left)){
        sf::Vector2i pos = sf::Mouse::getPosition(window);
        for (Slider * s : sliders){
          s->drag(pos.x,resY-pos.y);
        }
      }

      if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Left){
        for (Slider * s : sliders){
          s->mouseUp();
        }
      }

      // arrow keys nudge whichever slider was last clicked/dragged
      if (event.type == sf::Event::KeyPressed && activeSlider != nullptr){
        if (event.key.code == sf::Keyboard::Left || event.key.code == sf::Keyboard::Down){
          activeSlider->setPosition(activeSlider->getPosition()-0.01f);
        }
        else if (event.key.code == sf::Keyboard::Right || event.key.code == sf::Keyboard::Up){
          activeSlider->setPosition(activeSlider->getPosition()+0.01f);
        }
      }

    }

    window.clear(sf::Color::White);
    glClearColor(1.0f,1.0f,1.0f,1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    physClock.restart();

    if (randomizeButton.getState()){
      particles.randomisePositions();
      randomizeButton.setState(false);
    }

    // now for some very inelegant slider logic!
    int n = std::floor(particlesSlider.getPosition()*float(N));
    particlesSlider.setLabel("Particles: "+fixedLengthNumber(n,4));
    int nPdNi = particles.size()-particles.getHCount();
    if (n < nPdNi){
      while (n < particles.size()-particles.getHCount()){
        particles.removePdNiParticle();
      }
    }
    else if (n > nPdNi){
      int i = nPdNi-1;
      while (i < n){
        particles.addParticle();
        i++;
      }
    }

    double val = proportionBigSlider.getPosition();
    proportionBigSlider.setLabel("Pd:Ni Ratio: "+fixedLengthNumber(val,4)+" Pd");

    if (val != propBig){
      particles.randomise(val);
      propBig = val;
    }

    // proteium is added on top of the Pd/Ni population: value is the
    //  ratio of H particles to Pd particles (0 to MAX_H_RATIO, 1 = 1:1)
    val = proteiumSlider.getPosition()*MAX_H_RATIO;
    proteiumSlider.setLabel("Proteium Ratio: "+fixedLengthNumber(val,4)+" : 1 Pd");

    int currentPdNi = particles.size()-particles.getHCount();
    int nPd = std::floor(propBig*currentPdNi);
    particles.setHCount(uint64_t(std::round(val*nPd)));

    val = 1.0+massRatioSlider.getPosition()*maxMassRatio;
    particles.setMassRatio(val);
    massRatioSlider.setLabel("Mass Ratio: "+fixedLengthNumber(val,4));

    val = 1.0+radiusRatioSlider.getPosition()*maxRadiusRatio;
    particles.setRadiusRatio(val);
    radiusRatioSlider.setLabel("Size Ratio: "+fixedLengthNumber(val,4));

    val = pdAttractionSlider.getPosition()*maxPdAttraction;
    particles.setPdAttraction(val);
    pdAttractionSlider.setLabel("Pd Attraction: "+fixedLengthNumber(val,4));

    val = std::max(0.005,double(shakerSlider.getPosition())*shakerMaxPeriod);
    particles.setShakerPeriod(val);
    shakerSlider.setLabel("Shaker Period: "+fixedLengthNumber(val,5));

    val = amplitudeSlider.getPosition()*maxAmplitude;
    particles.setShakerAmplitude(val);
    amplitudeSlider.setLabel("Shaker Amplitude: "+fixedLengthNumber(val,4));

    val = std::min(std::max(0.1,double(restitutionSlider.getPosition())),0.98);
    particles.setCoeffientOfRestitution(val);
    restitutionSlider.setLabel("Coef. Restitution: "+fixedLengthNumber(val,4));
    particles.setTimeStep(dt*speed);

    if (!pause){
      for (int s = 0; s < subSamples; s++){
        particles.step();
        pRender.updatedTrack(particles);
      }
    }

    physDeltas[frameId] = physClock.getElapsedTime().asSeconds();

    renderClock.restart();

    glm::mat4 proj = camera.getVP();

    // restrict the 3D scene to the left of the control column, and drop
    //  it down from the top of the window by simTopMargin
    glViewport(0,0,simWidth,simHeight);
    glEnable(GL_DEPTH_TEST);
    pRender.setProjection(proj);
    pRender.draw(
      particles,
      frameId,
      simWidth,
      simHeight
    );
    glDisable(GL_DEPTH_TEST);

    // UI overlay below is 2D screen-space, full window, and always on top
    glViewport(0,0,resX,resY);

    if (newRecording.getState()){
      record.newFile();
      isRecording = !isRecording;
      newRecording.setState(false);
    }

    record.setSpeed(speed);

    if (frameId % saveFrequency == 0 && isRecording) {record.takeReading(particles);}

    if (isRecording){
      textRenderer.renderText(
        OD,
        "Recording to file:\n    "+record.fileName(),
        resX-400.0,resY-64.0*12,
        0.25f,
        glm::vec3(0.0f,0.0f,0.0f)
      );
    }

    if (debug){
      double delta = 0.0;
      double renderDelta = 0.0;
      double physDelta = 0.0;
      for (int n = 0; n < 60; n++){
        delta += deltas[n];
        renderDelta += renderDeltas[n];
        physDelta += physDeltas[n];
      }
      delta /= 60.0;
      renderDelta /= 60.0;
      physDelta /= 60.0;
      std::stringstream debugText;

      sf::Vector2i mouse = sf::Mouse::getPosition(window);

      glm::vec3 eye = camera.getEye();

      debugText << "Particles: " << N <<
        "\n" <<
        "Delta: " << fixedLengthNumber(delta,6) <<
        " (FPS: " << fixedLengthNumber(1.0/delta,4) << ")" <<
        "\n" <<
        "Render/Physics: " << fixedLengthNumber(renderDelta,6) << "/" << fixedLengthNumber(physDelta,6) <<
        "\n" <<
        "Mouse (" << fixedLengthNumber(mouse.x,4) << "," << fixedLengthNumber(mouse.y,4) << ")" <<
        "\n" <<
        "Camera eye (" << fixedLengthNumber(eye.x,4) << ", " << fixedLengthNumber(eye.y,4) << ", " << fixedLengthNumber(eye.z,4) << ")" <<
        "\n" <<
        "Order: " << fixedLengthNumber(particles.orderParameter(),6) << "\n";
      textRenderer.renderText(
        OD,
        debugText.str(),
        64.0f,resY-64.0f,
        0.5f,
        glm::vec3(0.0f,0.0f,0.0f)
      );
    }

    // more inelegant slider drawing
    for (Slider * s : sliders){
      s->draw(textRenderer,OD);
    }

    popups.draw(
      textRenderer,
      OD,
      clock.getElapsedTime().asSeconds()
    );

    newRecording.draw(
      textRenderer,
      OD
    );

    randomizeButton.draw(
      textRenderer,
      OD
    );

    if(pause){
      textRenderer.renderText(
        OD,
        "Space to resume",
        resX/3.,
        resY/2.,
        0.5,
        glm::vec3(0.,0.,0.),
        1.0
      );
    }

    window.display();

    deltas[frameId] = clock.getElapsedTime().asSeconds();
    renderDeltas[frameId] = renderClock.getElapsedTime().asSeconds();

    clock.restart();

    if (frameId == 60){
      frameId = 0;
    }
    else{
      frameId++;
    }
  }

  return 0;
}
