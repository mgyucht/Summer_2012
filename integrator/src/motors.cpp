// motors.cpp
// ----------
//
// motors.cpp contains the methods relating to force-dipole motor-network
// interactions in the integrator code.

#include "motors.h"

void Motors::step_motors()
{
    for (int i = 0; i < netSize; i++)
    {
        for (int j = 0; j < netSize; j++)
        {
            for (int k = 1; k < 4; k++)
            {
                // Two cases. Either motortimes[i] is positive or negative. No entry
                // can be zero.

                if (motortimes[(i * netSize + j) * 3 + k - 1] > 0) // The motor is bound to the network
                {
                    // If this motor is about to become unbound from the network
                    if (motortimes[(i * netSize + j) * 3 + k - 1] <= TIMESTEP) 
                    {
                        motortimes[(i * netSize + j) * 3 + k - 1] = generate_unbound_time();
                    } else
                    {
                        motortimes[(i * netSize + j) * 3 + k - 1] -= TIMESTEP;
                    }
                } else // The motor is not bound to the network
                {
                    if (motortimes[(i * netSize + j) * 3 + k - 1] >= -TIMESTEP)
                    {
                        if (std::abs(spr[i][j][k - 1]) < 1e-15) 
                        {
                            motortimes[(i * netSize + j) * 3 + k - 1] = generate_unbound_time();
                        } else
                        {
                            motortimes[(i * netSize + j) * 3 + k - 1] = generate_bound_time();
                        }
                    } else
                    {
                        motortimes[(i * netSize + j) * 3 + k - 1] += TIMESTEP;
                    }
                }
            }
        }
    }
}

double Motors::generate_bound_time()
{
    double mean = 4.0;
    return -log(1 - randDouble(0, 1)) / mean;
}

double Motors::generate_unbound_time()
{
    double mean = 20.0;
    return -log(1 - randDouble(0, 1)) / mean;
}

double Motors::getforce(int i, int j, int k)
{
    return motortimes[(i * netSize + j) * 3 + k] > 0 ? MOTORFORCE : 0;
}
