#ifndef _NONAFFINITY_H_
#define _NONAFFINITY_H_

// nonaffinity.h
// -------------
//
// Header file for nonaffinity.cpp.
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

double nonAffinity(const char * argv);

#endif /* _NONAFFINITY_H_ */
