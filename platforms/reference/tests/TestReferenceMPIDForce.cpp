/* -------------------------------------------------------------------------- *
 *                                   OpenMMMPID                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2015 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,  *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/**
 * This tests the Reference implementation of ReferenceMPIDForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMMPID.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/Units.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/MPIDForce.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/Vec3.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#define ASSERT_EQUAL_TOL_MOD(expected, found, tol, testname) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC_MOD(expected, found, tol,testname) {ASSERT_EQUAL_TOL_MOD((expected)[0], (found)[0], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[1], (found)[1], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[2], (found)[2], (tol),(testname));};


using namespace OpenMM;
using namespace std;

extern "C" OPENMM_EXPORT void registerMPIDReferenceKernelFactories();

const double TOL = 1e-4;

void make_waterbox(int natoms, double boxEdgeLength, MPIDForce *forceField,  vector<Vec3> &positions,
                   vector<pair<int, int> >& bonds, System &system) {
    bool do_charge = true;
    bool do_dpole  = true;
    bool do_qpole  = true;
    bool do_opole  = true;
    bool do_pol    = true;
    std::map < std::string, double > tholemap;
    std::map < std::string, double > polarmap;
    std::map < std::string, double > chargemap;
    std::map < std::string, std::vector<double> > dipolemap;
    std::map < std::string, std::vector<double> > quadrupolemap;
    std::map < std::string, std::vector<double> > octopolemap;
    std::map < std::string, MPIDForce::MultipoleAxisTypes > axesmap;
    std::map < std::string, std::vector<int> > anchormap;
    std::map < std::string, double > massmap;
    std::map < std::string, std::vector<int> > polgrpmap;
    std::map < std::string, std::vector<int> > cov12map;
    std::map < std::string, std::vector<int> > cov13map;

    axesmap["O"]  = MPIDForce::Bisector;
    axesmap["H1"] = MPIDForce::ZThenX;
    axesmap["H2"] = MPIDForce::ZThenX;

    chargemap["O"]  = -0.51966;
    chargemap["H1"] = 0.25983;
    chargemap["H2"] = 0.25983;
    if(!do_charge){
        chargemap["O"]  = 0.0;
        chargemap["H1"] = 0.0;
        chargemap["H2"] = 0.0;
    }

    int oanc[3] = {1, 2, 0};
    int h1anc[3] = {-1, 1, 0};
    int h2anc[3] = {-2, -1, 0};
    std::vector<int> oancv(&oanc[0], &oanc[3]);
    std::vector<int> h1ancv(&h1anc[0], &h1anc[3]);
    std::vector<int> h2ancv(&h2anc[0], &h2anc[3]);
    anchormap["O"]  = oancv;
    anchormap["H1"] = h1ancv;
    anchormap["H2"] = h2ancv;

    double od[3] = {0.0, 0.0, 0.00755612136146};
    double hd[3] = {-0.00204209484795, 0.0, -0.00307875299958};
    std::vector<double> odv(&od[0], &od[3]);
    std::vector<double> hdv(&hd[0], &hd[3]);
    if(!do_dpole){
        odv.assign(3, 0);
        hdv.assign(3, 0);
    }
    dipolemap["O"]  = odv;
    dipolemap["H1"] = hdv;
    dipolemap["H2"] = hdv;

    double oq[9] = {0.000354030721139,  0.0, 0.0,
                    0.0, -0.000390257077096, 0.0,
                    0.0, 0.0,  3.62263559571e-05};
    double hq[9] = {-3.42848248983e-05, 0.0, -1.89485963908e-06,
                     0.0,          -0.000100240875193,      0.0,
                    -1.89485963908e-06, 0.0,  0.000134525700091};

    std::vector<double> oqv(&oq[0], &oq[9]);
    std::vector<double> hqv(&hq[0], &hq[9]);
    if(!do_qpole){
        oqv.assign(9, 0);
        hqv.assign(9, 0);
    }
    quadrupolemap["O"]  = oqv;
    quadrupolemap["H1"] = hqv;
    quadrupolemap["H2"] = hqv;

    double oo[27] = {
                           0,                    0, -6.285758282686837e-07,
                           0,                    0,                    0,
        -6.285758282686837e-07,                    0,                    0,
                           0,                    0,                    0,
                           0,                    0, -9.452653225954594e-08,
                           0, -9.452653225954594e-08,                    0,
        -6.285758282686837e-07,                    0,                    0,
                           0, -9.452653225954594e-08,                    0,
                           0,                    0, 7.231018665791977e-07
    };
    double ho[27] = {
        -2.405600937552608e-07,                    0, -1.152422607026746e-06,
                           0, -6.415084018183151e-08,                    0,
        -1.152422607026746e-06,                    0, 3.047102424084479e-07,
                           0, -6.415084018183151e-08,                    0,
        -6.415084018183151e-08,                    0, -2.558537436767218e-06,
                           0, -2.558537436767218e-06,                    0,
        -1.152422607026746e-06,                    0, 3.047102424084479e-07,
                           0, -2.558537436767218e-06,                    0,
        3.047102424084479e-07,                    0, 3.710960043793964e-06
    };
    std::vector<double> oov(&oo[0], &oo[27]);
    std::vector<double> hov(&ho[0], &ho[27]);
    if(!do_opole){
        oov.assign(27, 0);
        hov.assign(27, 0);
    }
    octopolemap["O"]  = oov;
    octopolemap["H1"] = hov;
    octopolemap["H2"] = hov;

    polarmap["O"]  = 0.3069876538;
    polarmap["H1"] = 0.2813500172;
    polarmap["H2"] = 0.2813500172;

    polarmap["O"]  = 0.000837;
    polarmap["H1"] = 0.000496;
    polarmap["H2"] = 0.000496;

    tholemap["O"]  = 0.3900;
    tholemap["H1"] = 0.3900;
    tholemap["H2"] = 0.3900;

    massmap["O"]  = 15.999;
    massmap["H1"] = 1.0080000;
    massmap["H2"] = 1.0080000;

    int opg[3] = {0,1,2};
    int h1pg[3] = {-1,0,1};
    int h2pg[3] = {-2,-1,0};
    std::vector<int> opgv(&opg[0], &opg[3]);
    std::vector<int> h1pgv(&h1pg[0], &h1pg[3]);
    std::vector<int> h2pgv(&h2pg[0], &h2pg[3]);
    polgrpmap["O"] = opgv;
    polgrpmap["H1"] = h1pgv;
    polgrpmap["H2"] = h2pgv;

    int cov12o[2] = {1,2};
    int cov12h1[1] = {-1};
    int cov12h2[1] = {-2};
    std::vector<int> cov12ov(&cov12o[0], &cov12o[2]);
    std::vector<int> cov12h1v(&cov12h1[0], &cov12h1[1]);
    std::vector<int> cov12h2v(&cov12h2[0], &cov12h2[1]);
    cov12map["O"] = cov12ov;
    cov12map["H1"] = cov12h1v;
    cov12map["H2"] = cov12h2v;

    int cov13h1[1] = {1};
    int cov13h2[1] = {-1};
    std::vector<int> cov13h1v(&cov13h1[0], &cov13h1[1]);
    std::vector<int> cov13h2v(&cov13h2[0], &cov13h2[1]);
    cov13map["O"] = std::vector<int>();
    cov13map["H1"] = cov13h1v;
    cov13map["H2"] = cov13h2v;
    positions.clear();
    if (natoms == 6) {
        const double coords[6][3] = {
            {  2.000000, 2.000000, 2.000000},
            {  2.500000, 2.000000, 3.000000},
            {  1.500000, 2.000000, 3.000000},
            {  0.000000, 0.000000, 0.000000},
            {  0.500000, 0.000000, 1.000000},
            { -0.500000, 0.000000, 1.000000}
        };
        for (int atom = 0; atom < natoms; ++atom)
            positions.push_back(Vec3(coords[atom][0], coords[atom][1], coords[atom][2])*OpenMM::NmPerAngstrom);
    }
    else if (natoms == 375) {
        const double coords[375][3] = {
            { -6.22, -6.25, -6.24 },
            { -5.32, -6.03, -6.00 },
            { -6.75, -5.56, -5.84 },
            { -3.04, -6.23, -6.19 },
            { -3.52, -5.55, -5.71 },
            { -3.59, -6.43, -6.94 },
            {  0.02, -6.23, -6.14 },
            { -0.87, -5.97, -6.37 },
            {  0.53, -6.03, -6.93 },
            {  3.10, -6.20, -6.27 },
            {  3.87, -6.35, -5.72 },
            {  2.37, -6.11, -5.64 },
            {  6.18, -6.14, -6.20 },
            {  6.46, -6.66, -5.44 },
            {  6.26, -6.74, -6.94 },
            { -6.21, -3.15, -6.24 },
            { -6.23, -3.07, -5.28 },
            { -6.02, -2.26, -6.55 },
            { -3.14, -3.07, -6.16 },
            { -3.38, -3.63, -6.90 },
            { -2.18, -3.05, -6.17 },
            { -0.00, -3.16, -6.23 },
            { -0.03, -2.30, -6.67 },
            {  0.05, -2.95, -5.29 },
            {  3.08, -3.11, -6.14 },
            {  2.65, -2.55, -6.79 },
            {  3.80, -3.53, -6.62 },
            {  6.16, -3.14, -6.16 },
            {  7.04, -3.32, -6.51 },
            {  5.95, -2.27, -6.51 },
            { -6.20, -0.04, -6.15 },
            { -5.43,  0.32, -6.59 },
            { -6.95,  0.33, -6.62 },
            { -3.10, -0.06, -6.19 },
            { -3.75,  0.42, -6.69 },
            { -2.46,  0.60, -5.93 },
            {  0.05, -0.01, -6.17 },
            { -0.10,  0.02, -7.12 },
            { -0.79,  0.16, -5.77 },
            {  3.03,  0.00, -6.19 },
            {  3.54,  0.08, -7.01 },
            {  3.69, -0.22, -5.53 },
            {  6.17,  0.05, -6.19 },
            {  5.78, -0.73, -6.57 },
            {  7.09, -0.17, -6.04 },
            { -6.20,  3.15, -6.25 },
            { -6.59,  3.18, -5.37 },
            { -5.87,  2.25, -6.33 },
            { -3.09,  3.04, -6.17 },
            { -3.88,  3.58, -6.26 },
            { -2.41,  3.54, -6.63 },
            {  0.00,  3.06, -6.26 },
            { -0.71,  3.64, -6.00 },
            {  0.65,  3.15, -5.55 },
            {  3.14,  3.06, -6.23 },
            {  3.11,  3.31, -5.30 },
            {  2.38,  3.49, -6.63 },
            {  6.19,  3.14, -6.25 },
            {  6.82,  3.25, -5.54 },
            {  5.76,  2.30, -6.07 },
            { -6.22,  6.26, -6.19 },
            { -6.22,  5.74, -7.00 },
            { -5.89,  5.67, -5.52 },
            { -3.04,  6.24, -6.20 },
            { -3.08,  5.28, -6.17 },
            { -3.96,  6.52, -6.25 },
            { -0.05,  6.21, -6.16 },
            {  0.82,  6.58, -6.06 },
            {  0.01,  5.64, -6.93 },
            {  3.10,  6.25, -6.15 },
            {  3.64,  5.47, -6.31 },
            {  2.46,  6.24, -6.87 },
            {  6.22,  6.20, -6.27 },
            {  5.37,  6.42, -5.88 },
            {  6.80,  6.07, -5.51 },
            { -6.19, -6.15, -3.13 },
            { -6.37, -7.01, -3.51 },
            { -6.25, -6.29, -2.18 },
            { -3.10, -6.27, -3.11 },
            { -2.29, -5.77, -2.99 },
            { -3.80, -5.62, -2.98 },
            { -0.03, -6.18, -3.15 },
            { -0.07, -7.05, -2.75 },
            {  0.68, -5.74, -2.70 },
            {  3.10, -6.14, -3.07 },
            {  2.35, -6.72, -3.23 },
            {  3.86, -6.65, -3.37 },
            {  6.22, -6.20, -3.16 },
            {  6.82, -6.36, -2.43 },
            {  5.35, -6.13, -2.75 },
            { -6.26, -3.13, -3.12 },
            { -6.16, -2.27, -2.70 },
            { -5.36, -3.47, -3.18 },
            { -3.11, -3.05, -3.14 },
            { -3.31, -3.96, -3.34 },
            { -2.77, -3.06, -2.24 },
            {  0.00, -3.13, -3.16 },
            {  0.48, -2.37, -2.81 },
            { -0.57, -3.40, -2.44 },
            {  3.09, -3.09, -3.16 },
            {  2.41, -3.19, -2.49 },
            {  3.91, -3.07, -2.67 },
            {  6.19, -3.04, -3.08 },
            {  5.64, -3.61, -3.61 },
            {  6.93, -3.58, -2.82 },
            { -6.18, -0.00, -3.04 },
            { -6.00, -0.59, -3.78 },
            { -6.79,  0.64, -3.39 },
            { -3.05, -0.03, -3.07 },
            { -2.95,  0.80, -3.52 },
            { -4.00, -0.20, -3.07 },
            { -0.03,  0.03, -3.06 },
            { -0.33, -0.37, -3.87 },
            {  0.89, -0.21, -2.99 },
            {  3.13, -0.05, -3.10 },
            {  3.44,  0.81, -3.34 },
            {  2.21,  0.07, -2.86 },
            {  6.20, -0.05, -3.13 },
            {  6.89,  0.60, -3.20 },
            {  5.58,  0.30, -2.49 },
            { -6.23,  3.09, -3.16 },
            { -5.62,  3.79, -2.94 },
            { -6.33,  2.60, -2.33 },
            { -3.10,  3.08, -3.04 },
            { -3.84,  3.47, -3.51 },
            { -2.40,  3.01, -3.69 },
            {  0.01,  3.04, -3.11 },
            { -0.56,  3.59, -3.64 },
            {  0.28,  3.60, -2.38 },
            {  3.04,  3.11, -3.09 },
            {  3.49,  2.30, -2.87 },
            {  3.70,  3.66, -3.51 },
            {  6.15,  3.14, -3.11 },
            {  6.52,  2.52, -3.74 },
            {  6.72,  3.06, -2.34 },
            { -6.22,  6.15, -3.13 },
            { -5.49,  6.21, -2.51 },
            { -6.56,  7.04, -3.18 },
            { -3.11,  6.24, -3.05 },
            { -3.76,  5.83, -3.62 },
            { -2.26,  5.92, -3.37 },
            {  0.03,  6.25, -3.07 },
            {  0.34,  5.63, -3.73 },
            { -0.87,  6.00, -2.91 },
            {  3.07,  6.15, -3.08 },
            {  3.29,  6.92, -2.56 },
            {  3.39,  6.35, -3.96 },
            {  6.22,  6.14, -3.12 },
            {  5.79,  6.38, -2.29 },
            {  6.25,  6.96, -3.62 },
            { -6.21, -6.20, -0.06 },
            { -5.79, -6.87,  0.48 },
            { -6.43, -5.50,  0.54 },
            { -3.16, -6.21, -0.02 },
            { -2.50, -6.87,  0.20 },
            { -2.77, -5.37,  0.23 },
            { -0.00, -6.14, -0.00 },
            {  0.68, -6.72, -0.33 },
            { -0.64, -6.73,  0.38 },
            {  3.03, -6.20, -0.01 },
            {  3.77, -6.56, -0.51 },
            {  3.43, -5.85,  0.78 },
            {  6.25, -6.16, -0.00 },
            {  5.36, -6.09, -0.36 },
            {  6.24, -6.97,  0.49 },
            { -6.24, -3.05, -0.01 },
            { -6.35, -3.64,  0.73 },
            { -5.42, -3.33, -0.42 },
            { -3.09, -3.06,  0.05 },
            { -2.44, -3.62, -0.38 },
            { -3.90, -3.21, -0.43 },
            {  0.05, -3.10,  0.02 },
            { -0.31, -2.35, -0.43 },
            { -0.63, -3.77,  0.01 },
            {  3.05, -3.09, -0.04 },
            {  3.28, -3.90,  0.41 },
            {  3.65, -2.43,  0.30 },
            {  6.20, -3.04, -0.03 },
            {  5.66, -3.31,  0.71 },
            {  6.78, -3.79, -0.19 },
            { -6.18,  0.04, -0.04 },
            { -6.73, -0.73, -0.15 },
            { -5.98,  0.06,  0.89 },
            { -3.11, -0.04, -0.04 },
            { -3.36, -0.08,  0.87 },
            { -2.70,  0.81, -0.14 },
            { -0.02, -0.02, -0.05 },
            { -0.45,  0.28,  0.75 },
            {  0.90,  0.15,  0.07 },
            {  3.04,  0.02, -0.01 },
            {  3.26, -0.82,  0.38 },
            {  3.89,  0.45, -0.13 },
            {  6.19,  0.05, -0.03 },
            {  5.52, -0.56,  0.25 },
            {  7.01, -0.29,  0.32 },
            { -6.14,  3.08,  0.00 },
            { -6.83,  2.82,  0.61 },
            { -6.59,  3.64, -0.64 },
            { -3.05,  3.09, -0.04 },
            { -3.79,  2.50,  0.09 },
            { -3.18,  3.80,  0.59 },
            {  0.02,  3.14,  0.04 },
            { -0.89,  3.04, -0.19 },
            {  0.49,  2.57, -0.57 },
            {  3.14,  3.15,  0.00 },
            {  3.28,  2.28,  0.37 },
            {  2.30,  3.08, -0.45 },
            {  6.27,  3.08, -0.00 },
            {  5.55,  2.54, -0.33 },
            {  5.83,  3.87,  0.34 },
            { -6.18,  6.15, -0.03 },
            { -6.45,  6.21,  0.88 },
            { -6.26,  7.05, -0.36 },
            { -3.06,  6.19, -0.05 },
            { -2.84,  6.64,  0.76 },
            { -3.99,  5.96,  0.03 },
            { -0.00,  6.20,  0.06 },
            { -0.67,  5.99, -0.59 },
            {  0.76,  6.46, -0.44 },
            {  3.10,  6.26, -0.03 },
            {  3.57,  6.09,  0.78 },
            {  2.57,  5.47, -0.18 },
            {  6.26,  6.18,  0.02 },
            {  5.53,  5.64, -0.29 },
            {  5.95,  7.08, -0.06 },
            { -6.26, -6.21,  3.07 },
            { -5.98, -6.38,  3.97 },
            { -5.46, -5.94,  2.62 },
            { -3.10, -6.24,  3.04 },
            { -2.69, -6.51,  3.87 },
            { -3.43, -5.35,  3.21 },
            { -0.03, -6.16,  3.06 },
            {  0.83, -6.00,  3.42 },
            { -0.30, -6.99,  3.45 },
            {  3.15, -6.25,  3.11 },
            {  2.77, -5.60,  3.72 },
            {  2.68, -6.10,  2.28 },
            {  6.20, -6.21,  3.16 },
            {  5.75, -6.73,  2.50 },
            {  6.69, -5.56,  2.66 },
            { -6.17, -3.10,  3.04 },
            { -6.82, -2.44,  3.28 },
            { -6.12, -3.69,  3.80 },
            { -3.08, -3.04,  3.11 },
            { -3.59, -3.56,  3.72 },
            { -2.97, -3.61,  2.34 },
            {  0.01, -3.04,  3.11 },
            { -0.86, -3.41,  3.20 },
            {  0.56, -3.78,  2.86 },
            {  3.07, -3.07,  3.15 },
            {  3.81, -3.68,  3.13 },
            {  2.80, -2.98,  2.23 },
            {  6.20, -3.04,  3.13 },
            {  5.48, -3.64,  2.92 },
            {  6.98, -3.49,  2.81 },
            { -6.18, -0.05,  3.12 },
            { -6.41,  0.66,  3.69 },
            { -6.33,  0.28,  2.23 },
            { -3.05,  0.03,  3.10 },
            { -3.46, -0.42,  3.83 },
            { -3.57, -0.19,  2.33 },
            {  0.03, -0.02,  3.15 },
            {  0.23, -0.08,  2.21 },
            { -0.81,  0.41,  3.18 },
            {  3.09,  0.00,  3.03 },
            {  2.48, -0.29,  3.71 },
            {  3.91,  0.16,  3.51 },
            {  6.19, -0.06,  3.11 },
            {  6.05,  0.47,  2.33 },
            {  6.59,  0.52,  3.74 },
            { -6.20,  3.05,  3.05 },
            { -6.87,  3.73,  3.17 },
            { -5.55,  3.24,  3.73 },
            { -3.11,  3.06,  3.15 },
            { -3.64,  3.74,  2.71 },
            { -2.32,  3.00,  2.62 },
            {  0.02,  3.05,  3.06 },
            { -0.87,  3.14,  3.38 },
            {  0.48,  3.82,  3.42 },
            {  3.07,  3.10,  3.16 },
            {  3.95,  3.44,  2.97 },
            {  2.76,  2.73,  2.32 },
            {  6.19,  3.07,  3.16 },
            {  7.02,  3.30,  2.72 },
            {  5.52,  3.27,  2.51 },
            { -6.19,  6.24,  3.15 },
            { -5.56,  5.88,  2.52 },
            { -7.05,  5.96,  2.83 },
            { -3.10,  6.14,  3.08 },
            { -2.34,  6.69,  3.27 },
            { -3.86,  6.69,  3.29 },
            { -0.04,  6.24,  3.13 },
            {  0.63,  6.54,  2.53 },
            {  0.08,  5.29,  3.18 },
            {  3.12,  6.24,  3.14 },
            {  3.57,  5.82,  2.40 },
            {  2.23,  5.90,  3.12 },
            {  6.25,  6.19,  3.06 },
            {  5.55,  5.59,  3.32 },
            {  6.08,  6.99,  3.55 },
            { -6.20, -6.16,  6.15 },
            { -6.29, -5.99,  7.09 },
            { -6.09, -7.11,  6.09 },
            { -3.09, -6.19,  6.27 },
            { -2.56, -5.90,  5.52 },
            { -3.80, -6.69,  5.87 },
            {  0.02, -6.25,  6.24 },
            { -0.70, -5.70,  6.51 },
            {  0.25, -5.93,  5.36 },
            {  3.11, -6.18,  6.14 },
            {  3.76, -6.54,  6.74 },
            {  2.29, -6.20,  6.64 },
            {  6.22, -6.17,  6.15 },
            {  6.61, -6.98,  6.47 },
            {  5.56, -5.94,  6.81 },
            { -6.21, -3.10,  6.14 },
            { -6.76, -2.66,  6.78 },
            { -5.51, -3.50,  6.65 },
            { -3.13, -3.05,  6.18 },
            { -2.19, -3.14,  6.34 },
            { -3.50, -3.89,  6.43 },
            {  0.01, -3.06,  6.15 },
            { -0.06, -2.81,  7.07 },
            { -0.25, -3.98,  6.13 },
            {  3.04, -3.09,  6.17 },
            {  3.84, -3.51,  5.84 },
            {  3.25, -2.85,  7.08 },
            {  6.26, -3.13,  6.19 },
            {  6.01, -2.20,  6.09 },
            {  5.47, -3.55,  6.54 },
            { -6.20,  0.01,  6.27 },
            { -5.79, -0.70,  5.78 },
            { -6.67,  0.51,  5.60 },
            { -3.13,  0.01,  6.14 },
            { -3.53, -0.35,  6.94 },
            { -2.21,  0.17,  6.39 },
            { -0.04, -0.04,  6.20 },
            {  0.26,  0.47,  5.46 },
            {  0.51,  0.22,  6.93 },
            {  3.10, -0.05,  6.23 },
            {  2.33,  0.44,  5.95 },
            {  3.85,  0.45,  5.92 },
            {  6.19, -0.01,  6.26 },
            {  7.05,  0.16,  5.88 },
            {  5.58,  0.02,  5.52 },
            { -6.22,  3.04,  6.17 },
            { -5.45,  3.57,  5.95 },
            { -6.62,  3.50,  6.92 },
            { -3.09,  3.16,  6.21 },
            { -3.71,  2.75,  5.61 },
            { -2.60,  2.43,  6.59 },
            { -0.02,  3.10,  6.26 },
            {  0.89,  3.27,  6.05 },
            { -0.44,  2.94,  5.41 },
            {  3.12,  3.04,  6.23 },
            {  2.31,  3.53,  6.43 },
            {  3.59,  3.60,  5.60 },
            {  6.23,  3.05,  6.24 },
            {  5.92,  3.91,  6.54 },
            {  6.02,  3.03,  5.30 },
            { -6.15,  6.21,  6.24 },
            { -6.27,  6.46,  5.32 },
            { -7.00,  5.85,  6.51 },
            { -3.07,  6.15,  6.22 },
            { -3.98,  6.27,  5.94 },
            { -2.66,  7.01,  6.10 },
            {  0.04,  6.20,  6.25 },
            { -0.38,  5.50,  5.75 },
            { -0.36,  7.00,  5.93 },
            {  3.12,  6.15,  6.24 },
            {  3.66,  6.88,  5.93 },
            {  2.25,  6.33,  5.86 },
            {  6.20,  6.27,  6.19 },
            {  5.46,  5.65,  6.19 },
            {  6.97,  5.73,  6.39 }
        };
        for (int atom = 0; atom < natoms; ++atom)
            positions.push_back(Vec3(coords[atom][0], coords[atom][1], coords[atom][2])*OpenMM::NmPerAngstrom);
    }
    else
        throw exception();

    system.setDefaultPeriodicBoxVectors(Vec3(boxEdgeLength, 0, 0),
                                        Vec3(0, boxEdgeLength, 0),
                                        Vec3(0, 0, boxEdgeLength));

    const char* atom_types[3] = {"O", "H1", "H2"};
    for(int atom = 0; atom < natoms; ++atom){
        const char* element = atom_types[atom%3];
        double damp = polarmap[element];
        if(!do_pol) damp = 0.0;
        double alpha = pow(damp, 1.0/6.0);
        int atomz = atom + anchormap[element][0];
        int atomx = atom + anchormap[element][1];
        int atomy = anchormap[element][2]==0 ? -1 : atom + anchormap[element][2];
        forceField->addMultipole(chargemap[element], dipolemap[element], quadrupolemap[element], octopolemap[element],
                                 axesmap[element], atomz, atomx, atomy, tholemap[element], alpha, damp);
        system.addParticle(massmap[element]);
        // Polarization groups
        std::vector<int> tmppol;
        std::vector<int>& polgrps = polgrpmap[element];
        for(int i=0; i < polgrps.size(); ++i)
            tmppol.push_back(polgrps[i]+atom);
        if(!tmppol.empty())
           forceField->setCovalentMap(atom, MPIDForce::PolarizationCovalent11, tmppol);
        // 1-2 covalent groups
        std::vector<int> tmp12;
        std::vector<int>& cov12s = cov12map[element];
        for(int i=0; i < cov12s.size(); ++i)
            tmp12.push_back(cov12s[i]+atom);
        if(!tmp12.empty())
           forceField->setCovalentMap(atom, MPIDForce::Covalent12, tmp12);
        // 1-3 covalent groups
        std::vector<int> tmp13;
        std::vector<int>& cov13s = cov13map[element];
        for(int i=0; i < cov13s.size(); ++i)
            tmp13.push_back(cov13s[i]+atom);
        if(!tmp13.empty())
           forceField->setCovalentMap(atom, MPIDForce::Covalent13, tmp13);
    }
}

static void check_finite_differences(vector<Vec3> analytic_forces, Context &context, vector<Vec3> positions)
{
    // Take a small step in the direction of the energy gradient and see whether the potential energy changes by the expected amount.

    double norm = 0.0;
    for (auto& f : analytic_forces)
        norm += f.dot(f);
    norm = std::sqrt(norm);
    const double stepSize = 1e-3;
    double step = 0.5*stepSize/norm;
    vector<Vec3> positions2(analytic_forces.size()), positions3(analytic_forces.size());
    for (int i = 0; i < (int) positions.size(); ++i) {
        Vec3 p = positions[i];
        Vec3 f = analytic_forces[i];
        positions2[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
        positions3[i] = Vec3(p[0]+f[0]*step, p[1]+f[1]*step, p[2]+f[2]*step);
    }
    context.setPositions(positions2);
    State state2 = context.getState(State::Energy);
    context.setPositions(positions3);
    State state3 = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state3.getPotentialEnergy())/stepSize, 1e-4);
}

void testWaterDimerEnergyAndForces() {
    const double cutoff = 6.0*OpenMM::NmPerAngstrom;
    double boxEdgeLength = 20*OpenMM::NmPerAngstrom;
    const double alpha = 3.0;
    const int grid = 64;
    MPIDForce* forceField = new MPIDForce();

    vector<Vec3> positions;
    vector<pair<int, int> > bonds;
    System system;

    const int numAtoms = 6;

    make_waterbox(numAtoms, boxEdgeLength, forceField,  positions, bonds, system);
    forceField->setNonbondedMethod(OpenMM::MPIDForce::PME);
    forceField->setPolarizationType(MPIDForce::Extrapolated);
    forceField->setPMEParameters(alpha, grid, grid, grid);
    forceField->setCutoffDistance(cutoff);
    system.addForce(forceField);

    VerletIntegrator integrator(0.01);
    Context context(system, integrator, Platform::getPlatformByName("Reference"));
    context.setPositions(positions);
    double refenergy =-2.298363777;
    vector<Vec3> refforces(6);
    refforces[0] = Vec3(  -137.7203283, -183.9079957,  32.80100273);
    refforces[1] = Vec3(   37.56519756, -4.813998574,  7.691308969);
    refforces[2] = Vec3(   40.99501057,  119.8240071,  34.06469001);
    refforces[3] = Vec3(  -117.5860009,  -101.736468, -32.73843882);
    refforces[4] = Vec3(   129.5884295,  167.1137166, -18.58460594);
    refforces[5] = Vec3(   47.15774644,  3.520730293, -23.23411526);

    State state = context.getState(State::Forces | State::Energy);
    double energy = state.getPotentialEnergy();
    const vector<Vec3>& forces = state.getForces();

    ASSERT_EQUAL_TOL(refenergy, energy, 1E-4);
    for (int n = 0; n < numAtoms; ++n)
        ASSERT_EQUAL_VEC(refforces[n], forces[n], 1E-4);

    // Just to be sure, check the finite difference gradients
    check_finite_differences(forces, context, positions);
}



int main(int numberOfArguments, char* argv[]) {

    try {
        std::cout << "TestReferenceMPIDForce running test..." << std::endl;
        registerMPIDReferenceKernelFactories();

        testWaterDimerEnergyAndForces();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
