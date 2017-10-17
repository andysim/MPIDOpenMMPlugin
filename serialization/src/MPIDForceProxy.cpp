/* -------------------------------------------------------------------------- *
 *                                OpenMMMPID                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2016 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
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

#include "openmm/serialization/MPIDForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/MPIDForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

MPIDForceProxy::MPIDForceProxy() : SerializationProxy("MPIDForce") {
}

static void getCovalentTypes(std::vector<std::string>& covalentTypes) {

    covalentTypes.push_back("Covalent12");
    covalentTypes.push_back("Covalent13");
    covalentTypes.push_back("Covalent14");
    covalentTypes.push_back("Covalent15");

    covalentTypes.push_back("PolarizationCovalent11");
    covalentTypes.push_back("PolarizationCovalent12");
    covalentTypes.push_back("PolarizationCovalent13");
    covalentTypes.push_back("PolarizationCovalent14");
}

static void addCovalentMap(SerializationNode& particleExclusions, int particleIndex, std::string mapName, std::vector< int > covalentMap) {
    SerializationNode& map   = particleExclusions.createChildNode(mapName);
    for (unsigned int ii = 0; ii < covalentMap.size(); ii++) {
        map.createChildNode("Cv").setIntProperty("v", covalentMap[ii]);
    }
}

void loadCovalentMap(const SerializationNode& map, std::vector< int >& covalentMap) {
    for (unsigned int ii = 0; ii < map.getChildren().size(); ii++) {
        covalentMap.push_back(map.getChildren()[ii].getIntProperty("v"));
    }
}

void MPIDForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 5);
    const MPIDForce& force = *reinterpret_cast<const MPIDForce*>(object);

    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setIntProperty("nonbondedMethod",                  force.getNonbondedMethod());
    node.setIntProperty("polarizationType",                 force.getPolarizationType());
    node.setIntProperty("mutualInducedMaxIterations",       force.getMutualInducedMaxIterations());

    node.setDoubleProperty("cutoffDistance",                force.getCutoffDistance());
    double alpha;
    int nx, ny, nz;
    force.getPMEParameters(alpha, nx, ny, nz);
    node.setDoubleProperty("aEwald",                        alpha);
    node.setDoubleProperty("mutualInducedTargetEpsilon",    force.getMutualInducedTargetEpsilon());
    node.setDoubleProperty("ewaldErrorTolerance",           force.getEwaldErrorTolerance());

    SerializationNode& gridDimensionsNode  = node.createChildNode("MultipoleParticleGridDimension");
    gridDimensionsNode.setIntProperty("d0", nx).setIntProperty("d1", ny).setIntProperty("d2", nz); 
    
    SerializationNode& coefficients = node.createChildNode("ExtrapolationCoefficients");
    vector<double> coeff = force.getExtrapolationCoefficients();
    for (int i = 0; i < coeff.size(); i++) {
        stringstream key;
        key << "c" << i;
        coefficients.setDoubleProperty(key.str(), coeff[i]);
    }

    std::vector<std::string> covalentTypes;
    getCovalentTypes(covalentTypes);

    SerializationNode& particles = node.createChildNode("MultipoleParticles");
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force.getNumMultipoles()); ii++) {

        int axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY;
        double charge, thole, dampingFactor, polarity;

        std::vector<double> molecularDipole;
        std::vector<double> molecularQuadrupole;
        std::vector<double> molecularOctopole;

        force.getMultipoleParameters(ii, charge, molecularDipole, molecularQuadrupole, molecularOctopole,
                                     axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);

        SerializationNode& particle    = particles.createChildNode("Particle");
        particle.setIntProperty("axisType", axisType).setIntProperty("multipoleAtomZ", multipoleAtomZ).setIntProperty("multipoleAtomX", multipoleAtomX).setIntProperty("multipoleAtomY", multipoleAtomY);
        particle.setDoubleProperty("charge", charge).setDoubleProperty("thole", thole).setDoubleProperty("damp", dampingFactor).setDoubleProperty("polarity", polarity);

        SerializationNode& dipole      = particle.createChildNode("Dipole");
        dipole.setDoubleProperty("d0", molecularDipole[0]).setDoubleProperty("d1", molecularDipole[1]).setDoubleProperty("d2", molecularDipole[2]);

        SerializationNode& quadrupole  = particle.createChildNode("Quadrupole");
        quadrupole.setDoubleProperty("q0", molecularQuadrupole[0]).setDoubleProperty("q1", molecularQuadrupole[1]).setDoubleProperty("q2", molecularQuadrupole[2]);
        quadrupole.setDoubleProperty("q3", molecularQuadrupole[3]).setDoubleProperty("q4", molecularQuadrupole[4]).setDoubleProperty("q5", molecularQuadrupole[5]);
        quadrupole.setDoubleProperty("q6", molecularQuadrupole[6]).setDoubleProperty("q7", molecularQuadrupole[7]).setDoubleProperty("q8", molecularQuadrupole[8]);

        SerializationNode& octopole  = particle.createChildNode("Octopole");
        octopole.setDoubleProperty("o0", molecularOctopole[0]).setDoubleProperty("o1", molecularOctopole[1]).setDoubleProperty("o2", molecularOctopole[2]);
        octopole.setDoubleProperty("o3", molecularOctopole[3]).setDoubleProperty("o4", molecularOctopole[4]).setDoubleProperty("o5", molecularOctopole[5]);
        octopole.setDoubleProperty("o6", molecularOctopole[6]).setDoubleProperty("o7", molecularOctopole[7]).setDoubleProperty("o8", molecularOctopole[8]);
        octopole.setDoubleProperty("o9", molecularOctopole[9]).setDoubleProperty("o10", molecularOctopole[10]).setDoubleProperty("o11", molecularOctopole[11]);
        octopole.setDoubleProperty("o12", molecularOctopole[12]).setDoubleProperty("o13", molecularOctopole[13]).setDoubleProperty("o14", molecularOctopole[14]);
        octopole.setDoubleProperty("o15", molecularOctopole[15]).setDoubleProperty("o16", molecularOctopole[16]).setDoubleProperty("o17", molecularOctopole[17]);
        octopole.setDoubleProperty("o18", molecularOctopole[18]).setDoubleProperty("o19", molecularOctopole[19]).setDoubleProperty("o20", molecularOctopole[20]);
        octopole.setDoubleProperty("o21", molecularOctopole[21]).setDoubleProperty("o22", molecularOctopole[22]).setDoubleProperty("o23", molecularOctopole[23]);
        octopole.setDoubleProperty("o24", molecularOctopole[24]).setDoubleProperty("o25", molecularOctopole[25]).setDoubleProperty("o26", molecularOctopole[26]);

        for (unsigned int jj = 0; jj < covalentTypes.size(); jj++) {
            std::vector< int > covalentMap;
            force.getCovalentMap(ii, static_cast<MPIDForce::CovalentType>(jj), covalentMap);
            addCovalentMap(particle, ii, covalentTypes[jj], covalentMap);
        }
    }
}

void* MPIDForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version < 0 || version > 5)
        throw OpenMMException("Unsupported version number");
    MPIDForce* force = new MPIDForce();

    try {
        if (version > 3)
            force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setNonbondedMethod(static_cast<MPIDForce::NonbondedMethod>(node.getIntProperty("nonbondedMethod")));
        if (version >= 2)
            force->setPolarizationType(static_cast<MPIDForce::PolarizationType>(node.getIntProperty("polarizationType")));
        force->setMutualInducedMaxIterations(node.getIntProperty("mutualInducedMaxIterations"));

        force->setCutoffDistance(node.getDoubleProperty("cutoffDistance"));
        force->setMutualInducedTargetEpsilon(node.getDoubleProperty("mutualInducedTargetEpsilon"));
        force->setEwaldErrorTolerance(node.getDoubleProperty("ewaldErrorTolerance"));

        const SerializationNode& gridDimensionsNode  = node.getChildNode("MultipoleParticleGridDimension");
        force->setPMEParameters(node.getDoubleProperty("aEwald"), gridDimensionsNode.getIntProperty("d0"), gridDimensionsNode.getIntProperty("d1"), gridDimensionsNode.getIntProperty("d2"));
    
        if (version >= 3) {
            const SerializationNode& coefficients = node.getChildNode("ExtrapolationCoefficients");
            vector<double> coeff;
            for (int i = 0; ; i++) {
                stringstream key;
                key << "c" << i;
                if (coefficients.getProperties().find(key.str()) == coefficients.getProperties().end())
                    break;
                coeff.push_back(coefficients.getDoubleProperty(key.str()));
            }
            force->setExtrapolationCoefficients(coeff);
        }
        std::vector<std::string> covalentTypes;
        getCovalentTypes(covalentTypes);

        const SerializationNode& particles = node.getChildNode("MultipoleParticles");
        for (unsigned int ii = 0; ii < particles.getChildren().size(); ii++) {

            const SerializationNode& particle = particles.getChildren()[ii];

            std::vector<double> molecularDipole;
            const SerializationNode& dipole = particle.getChildNode("Dipole");
            molecularDipole.push_back(dipole.getDoubleProperty("d0"));
            molecularDipole.push_back(dipole.getDoubleProperty("d1"));
            molecularDipole.push_back(dipole.getDoubleProperty("d2"));

            std::vector<double> molecularQuadrupole;
            const SerializationNode& quadrupole = particle.getChildNode("Quadrupole");
            molecularQuadrupole.push_back(quadrupole.getDoubleProperty("q0"));
            molecularQuadrupole.push_back(quadrupole.getDoubleProperty("q1"));
            molecularQuadrupole.push_back(quadrupole.getDoubleProperty("q2"));
            molecularQuadrupole.push_back(quadrupole.getDoubleProperty("q3"));
            molecularQuadrupole.push_back(quadrupole.getDoubleProperty("q4"));
            molecularQuadrupole.push_back(quadrupole.getDoubleProperty("q5"));
            molecularQuadrupole.push_back(quadrupole.getDoubleProperty("q6"));
            molecularQuadrupole.push_back(quadrupole.getDoubleProperty("q7"));
            molecularQuadrupole.push_back(quadrupole.getDoubleProperty("q8"));

            std::vector<double> molecularOctopole;
            if(version>=5) {
                const SerializationNode& octopole = particle.getChildNode("Octopole");
                molecularOctopole.push_back(octopole.getDoubleProperty("o0"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o1"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o2"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o3"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o4"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o5"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o6"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o7"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o8"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o9"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o10"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o11"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o12"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o13"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o14"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o15"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o16"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o17"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o18"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o19"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o20"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o21"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o22"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o23"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o24"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o25"));
                molecularOctopole.push_back(octopole.getDoubleProperty("o26"));
            } else {
                for(int i = 0; i < 27; ++i)
                    molecularOctopole.push_back(0.0);
            }
            force->addMultipole(particle.getDoubleProperty("charge"), molecularDipole, molecularQuadrupole, molecularOctopole,
                                particle.getIntProperty("axisType"),
                                particle.getIntProperty("multipoleAtomZ"),
                                particle.getIntProperty("multipoleAtomX"),
                                particle.getIntProperty("multipoleAtomY"),
                                particle.getDoubleProperty("thole"),
                                particle.getDoubleProperty("damp"), particle.getDoubleProperty("polarity"));

            // covalent maps 

            for (unsigned int jj = 0; jj < covalentTypes.size(); jj++) {
                std::vector< int > covalentMap;
                loadCovalentMap(particle.getChildNode(covalentTypes[jj]), covalentMap);
                force->setCovalentMap(ii, static_cast<MPIDForce::CovalentType>(jj), covalentMap);
            }
        }

    }
    catch (...) {
        delete force;
        throw;
    }

    return force;
}
