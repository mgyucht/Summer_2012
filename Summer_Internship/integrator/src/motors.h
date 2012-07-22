#ifndef MOTORS_H_
#define MOTORS_H_
// motors.h
// --------
//
// The header file corresponing to motors.c. motors.h provides the important
// methods for force-dipole motor-network interactions in the integrator
// simulation.

#include <cmath>
#include "utils.h"

extern int netSize;
extern double TIMESTEP;

const double MOTORFORCE = 1e-2; // WHAT?

class Motors
{
    public:
    double*** spr;
    
    Motors(double*** sspr) : spr(sspr) 
    {
        motortimes = new double[6 * netSize * netSize];
    }
    
    ~Motors() 
    {
    }
    
    // This moves the motors to the next time step.
    void step_motors();
    
    // These methods draw from a random distribution to determine how long a
    // given motor will stay attached to or removed from the network.
    double generate_bound_time();
    double generate_unbound_time();
    
    double getforce(int /* row */, int /* col */, int /* spr */);
    
    private:
    double* motortimes;

};

#endif /* MOTORS_H_ */
