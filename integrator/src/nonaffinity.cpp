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

double affxpos(int r, int c)
{
    return c + r / 2.0 + affdel / 2.0 * ((2.0 * r) / (netSize - 1.0) - 1.0);
}

double affypos(int r)
{
    return sqrt(3.0) / 2.0 * r;
}

double nonAffinity(double *position)
{
  double xval, yval, currentx, currenty, prefactor, sqrdisp = 0, nonaffinity = 0;

  if (std::abs(affdel) < 1E-15)
  {
    return 0.0;
  }

  prefactor = 1; // / ((float) netSize * netSize * strain * strain);

  for (int row = 0; row < netSize; row++) {
    for (int col = 0; col < netSize; col++) {
      currentx = position[(row * netSize + col) * 2];
      currenty = position[(row * netSize + col) * 2 + 1];

      xval = affxpos(row, col);
      yval = affypos(row);

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
  double dxval, sqrdisp = 0.0, nonaffinity = 0.0, currentdx, currentdy;
  double prefactor = 1; /* / ((float) netSize * netSize * str_rate * str_rate
      * TIMESTEP * TIMESTEP); */

  for (int row = 0; row < netSize; row++)
  {
    for (int col = 0; col < netSize; col++)
    {
      currentdx = delta[(row * netSize + col) * 2] / TIMESTEP;
      currentdy = delta[(row * netSize + col) * 2 + 1] / TIMESTEP;

      // dxval = 1 / 2.0 * ((2.0 * row) / (netSize - 1) - 1) * str_rate;
      dxval = affvx(row, str_rate);

      sqrdisp += (currentdx - dxval) * (currentdx - dxval) + currentdy * currentdy;
    }
  }

  nonaffinity = prefactor * sqrdisp;
  return nonaffinity;
}
