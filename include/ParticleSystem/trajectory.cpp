#include <ParticleSystem/trajectory.h>

void Trajectory::takeReading(ParticleSystem & p){
    std::vector<State> s(p.size(),State(0.0,0.0,0.0,0.0));
    for (int i = 0; i < p.size(); i++){
        s[i].x = p.state[i*3];
        s[i].y = p.state[i*3+1];
        s[i].z = p.state[i*3+2];
        s[i].radius = p.parameters[i*2];
        size += 1;
    }
    trajectory.push_back(s);
    shakerPosition.push_back(p.shakerDisplacement+p.shakerAmplitude);
    parameters.push_back(
        Parameters(
            p.massRatio,
            p.radiusRatio,
            p.shakerPeriod,
            p.shakerAmplitude,
            p.coefficientOfRestitution,
            p.size(),
            p.orderParameter(),
            speed
        )
    );

    if (size*(sizeof(State)) > CACHE_SIZE_BYTES){
        save();
        ioMode = std::ios_base::app;
        clear();
    }
}

void Trajectory::save(){

    #if WINDOWS
        // threading needs to be handled differently
        //   apparently
        threadedSave(trajectory);
    #else
        std::thread job(&Trajectory::threadedSave,this,this->trajectory);
        job.detach();
    #endif

}

void Trajectory::threadedSave(std::vector<std::vector<State>> trajectory){
    std::ofstream out(file,ioMode);
    std::ofstream outShaker("shaker-"+file,ioMode);
    std::ofstream outParam("parameters-"+file,ioMode);
    if (out.is_open() && outShaker.is_open() && outParam.is_open()){
        for (int t = 0; t < trajectory.size(); t++){
            for (int i = 0; i < trajectory[t].size(); i++){
                out << trajectory[t][i].x << ", "
                    << trajectory[t][i].y << ", "
                    << trajectory[t][i].z << ", "
                    << trajectory[t][i].radius
                    << "\n";
            }
            out << "\n";
            outShaker << shakerPosition[t] << ",";
            outParam << parameters[t].massRatio << ","
                     << parameters[t].radiusRatio << ","
                     << parameters[t].period << ","
                     << parameters[t].amplitude << ","
                     << parameters[t].restitution << ","
                     << parameters[t].nParticles << ","
                     << parameters[t].order << ","
                     << parameters[t].speed << "\n";

        }
        out.close(); 
        outShaker.close();
        outParam.close();
    }
    else{
        throw IOException("Cannot open file "+file);
    }
}