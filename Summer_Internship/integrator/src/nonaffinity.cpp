// nonaffinity.cpp
// ---------------
//
// nonaffinity.cpp takes as its input a position file and outputs the nonaffinity measure
// for the network. The nonaffinity measure of a network is a way of measuring how much
// a network deforming nonaffinely varies from the affine deformation. It is calculated
// by the formula:
//
//          1
//     Γ = --- Σ (u_i - uaff_i)²
//         Nγ²
//
// where Γ is the nonaffinity measure, N is the number of nodes in the network, γ is the
// strain on the network, u_i is the displacement of the node in the nonaffine
// deformation, and uaff_i the displacement ofthe same node in the affine one.

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>

#include "nonaffinity.h"

double nonAffinity(double *position)
{
  double xval, yval, currentx, currenty, prefactor, sqrdisp = 0,
         nonaffinity = 0;

  if (std::abs(strain) < 1E-15)
  {
    return 0.0;
  }

  prefactor = 1 / ((float) netSize * netSize * strain * strain);

  for (int row = 0; row < netSize; row++)
  {
    for (int col = 0; col < netSize; col++)
    {
      currentx = position[(row * netSize + col) * 2];
      currenty = position[(row * netSize + col) * 2 + 1];

      xval = col + row / 2.0 + sqrt(3) / 4 * strain * netSize
        * (2 * ((double) row - netSize) / (netSize + 1) + 1);
      yval = sqrt(3.0) / 2.0 * row;

      double temp = (currentx - xval) * (currentx - xval) + (currenty - yval)
        * (currenty - yval);
      sqrdisp += temp;
    }
  }

  nonaffinity = prefactor * sqrdisp;

  return nonaffinity;

}

double nonAffinity_dd(double *position, double *delta, double str_rate)
{
  // dyval is always 0; the affine prediction is that nodes don't move in the y
  // direction.
  double dxval, currentx, currenty, sqrdisp = 0.0, nonaffinity = 0.0;
  double currentdx, currentdy, lastx, lasty;
  double prefactor = 1 / ((float) netSize * netSize * str_rate * str_rate
      * TIMESTEP * TIMESTEP);

  for (int row = 0; row < netSize; row++)
  {
    for (int col = 0; col < netSize; col++)
    {
      currentx = position[(row * netSize + col) * 2];
      currenty = position[(row * netSize + col) * 2 + 1];

      currentdx = delta[(row * netSize + col) * 2];
      currentdy = delta[(row * netSize + col) * 2 + 1];

      lastx = currentx - currentdx;
      lasty = currenty - currentdy;

      dxval = str_rate * sqrt(3.0) / 4.0 * netSize * (1 - 2 * ((float) (netSize - row)) / (netSize + 1)) * TIMESTEP;

      sqrdisp += (currentdx - dxval) * (currentdx - dxval) + currentdy * currentdy;
    }
  }

  nonaffinity = prefactor * sqrdisp;
  return nonaffinity;
}
