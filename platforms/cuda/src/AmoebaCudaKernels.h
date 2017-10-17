#ifndef MPID_OPENMM_CUDAKERNELS_H_
#define MPID_OPENMM_CUDAKERNELS_H_

/* -------------------------------------------------------------------------- *
 *                              OpenMMMPID                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2015 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs, Peter Eastman                                    *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "openmm/mpidKernels.h"
#include "openmm/kernels.h"
#include "openmm/System.h"
#include "CudaArray.h"
#include "CudaContext.h"
#include "CudaSort.h"
#include <cufft.h>

namespace OpenMM {

class CudaCalcMPIDGeneralizedKirkwoodForceKernel;

/**
 * This kernel is invoked by MPIDBondForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMPIDBondForceKernel : public CalcMPIDBondForceKernel {
public:
    CudaCalcMPIDBondForceKernel(std::string name, 
                                          const Platform& platform,
                                          CudaContext& cu,
                                          const System& system);
    ~CudaCalcMPIDBondForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MPIDBondForce this kernel will be used for
     */
    void initialize(const System& system, const MPIDBondForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MPIDBondForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MPIDBondForce& force);
private:
    class ForceInfo;
    int numBonds;
    CudaContext& cu;
    const System& system;
    CudaArray* params;
};

/**
 * This kernel is invoked by MPIDAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMPIDAngleForceKernel : public CalcMPIDAngleForceKernel {
public:
    CudaCalcMPIDAngleForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMPIDAngleForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MPIDAngleForce this kernel will be used for
     */
    void initialize(const System& system, const MPIDAngleForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MPIDAngleForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MPIDAngleForce& force);
private:
    class ForceInfo;
    int numAngles;
    CudaContext& cu;
    const System& system;
    CudaArray* params;
};

/**
 * This kernel is invoked by MPIDInPlaneAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMPIDInPlaneAngleForceKernel : public CalcMPIDInPlaneAngleForceKernel {
public:
    CudaCalcMPIDInPlaneAngleForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMPIDInPlaneAngleForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MPIDInPlaneAngleForce this kernel will be used for
     */
    void initialize(const System& system, const MPIDInPlaneAngleForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MPIDInPlaneAngleForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MPIDInPlaneAngleForce& force);
private:
    class ForceInfo;
    int numAngles;
    CudaContext& cu;
    const System& system;
    CudaArray* params;
};

/**
 * This kernel is invoked by MPIDPiTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMPIDPiTorsionForceKernel : public CalcMPIDPiTorsionForceKernel {
public:
    CudaCalcMPIDPiTorsionForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMPIDPiTorsionForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MPIDPiTorsionForce this kernel will be used for
     */
    void initialize(const System& system, const MPIDPiTorsionForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MPIDPiTorsionForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MPIDPiTorsionForce& force);
private:
    class ForceInfo;
    int numPiTorsions;
    CudaContext& cu;
    const System& system;
    CudaArray* params;
};

/**
 * This kernel is invoked by MPIDStretchBendForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMPIDStretchBendForceKernel : public CalcMPIDStretchBendForceKernel {
public:
    CudaCalcMPIDStretchBendForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMPIDStretchBendForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MPIDStretchBendForce this kernel will be used for
     */
    void initialize(const System& system, const MPIDStretchBendForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MPIDStretchBendForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MPIDStretchBendForce& force);
private:
    class ForceInfo;
    int numStretchBends;
    CudaContext& cu;
    const System& system;
    CudaArray* params1; // Equilibrium values
    CudaArray* params2; // force constants
};

/**
 * This kernel is invoked by MPIDOutOfPlaneBendForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMPIDOutOfPlaneBendForceKernel : public CalcMPIDOutOfPlaneBendForceKernel {
public:
    CudaCalcMPIDOutOfPlaneBendForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMPIDOutOfPlaneBendForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MPIDOutOfPlaneBendForce this kernel will be used for
     */
    void initialize(const System& system, const MPIDOutOfPlaneBendForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MPIDOutOfPlaneBendForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MPIDOutOfPlaneBendForce& force);
private:
    class ForceInfo;
    int numOutOfPlaneBends;
    CudaContext& cu;
    const System& system;
    CudaArray* params;
};

/**
 * This kernel is invoked by MPIDTorsionTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMPIDTorsionTorsionForceKernel : public CalcMPIDTorsionTorsionForceKernel {
public:
    CudaCalcMPIDTorsionTorsionForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMPIDTorsionTorsionForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MPIDTorsionTorsionForce this kernel will be used for
     */
    void initialize(const System& system, const MPIDTorsionTorsionForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
private:
    class ForceInfo;
    int numTorsionTorsions;
    int numTorsionTorsionGrids;
    CudaContext& cu;
    const System& system;
    CudaArray* gridValues;
    CudaArray* gridParams;
    CudaArray* torsionParams;
};

/**
 * This kernel is invoked by MPIDMultipoleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMPIDMultipoleForceKernel : public CalcMPIDMultipoleForceKernel {
public:
    CudaCalcMPIDMultipoleForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMPIDMultipoleForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MPIDMultipoleForce this kernel will be used for
     */
    void initialize(const System& system, const MPIDMultipoleForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
     /**
     * Get the LabFrame dipole moments of all particles.
     * 
     * @param context    the Context for which to get the induced dipoles
     * @param dipoles    the induced dipole moment of particle i is stored into the i'th element
     */
    void getLabFramePermanentDipoles(ContextImpl& context, std::vector<Vec3>& dipoles);
    /**
     * Get the induced dipole moments of all particles.
     * 
     * @param context    the Context for which to get the induced dipoles
     * @param dipoles    the induced dipole moment of particle i is stored into the i'th element
     */
    void getInducedDipoles(ContextImpl& context, std::vector<Vec3>& dipoles);
    /**
     * Get the total dipole moments of all particles.
     * 
     * @param context    the Context for which to get the induced dipoles
     * @param dipoles    the induced dipole moment of particle i is stored into the i'th element
     */
    void getTotalDipoles(ContextImpl& context, std::vector<Vec3>& dipoles);
    /**
     * Execute the kernel to calculate the electrostatic potential
     *
     * @param context        the context in which to execute this kernel
     * @param inputGrid      input grid coordinates
     * @param outputElectrostaticPotential output potential 
     */
    void getElectrostaticPotential(ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                   std::vector< double >& outputElectrostaticPotential);

   /** 
     * Get the system multipole moments
     *
     * @param context      context
     * @param outputMultipoleMoments (charge,
     *                                dipole_x, dipole_y, dipole_z,
     *                                quadrupole_xx, quadrupole_xy, quadrupole_xz,
     *                                quadrupole_yx, quadrupole_yy, quadrupole_yz,
     *                                quadrupole_zx, quadrupole_zy, quadrupole_zz)
     */
    void getSystemMultipoleMoments(ContextImpl& context, std::vector<double>& outputMultipoleMoments);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MPIDMultipoleForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MPIDMultipoleForce& force);
    /**
     * Get the parameters being used for PME.
     * 
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    void getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;
private:
    class ForceInfo;
    class SortTrait : public CudaSort::SortTrait {
        int getDataSize() const {return 8;}
        int getKeySize() const {return 4;}
        const char* getDataType() const {return "int2";}
        const char* getKeyType() const {return "int";}
        const char* getMinKey() const {return "(-2147483647 - 1)";}
        const char* getMaxKey() const {return "2147483647";}
        const char* getMaxValue() const {return "make_int2(2147483647, 2147483647)";}
        const char* getSortKey() const {return "value.y";}
    };
    void initializeScaleFactors();
    void computeInducedField(void** recipBoxVectorPointer);
    bool iterateDipolesByDIIS(int iteration);
    void computeExtrapolatedDipoles(void** recipBoxVectorPointer);
    void ensureMultipolesValid(ContextImpl& context);
    template <class T, class T4, class M4> void computeSystemMultipoleMoments(ContextImpl& context, std::vector<double>& outputMultipoleMoments);
    int numMultipoles, maxInducedIterations, maxExtrapolationOrder;
    int fixedFieldThreads, inducedFieldThreads, electrostaticsThreads;
    int gridSizeX, gridSizeY, gridSizeZ;
    double alpha, inducedEpsilon;
    bool usePME, hasQuadrupoles, hasInitializedScaleFactors, hasInitializedFFT, multipolesAreValid, hasCreatedEvent;
    MPIDMultipoleForce::PolarizationType polarizationType;
    CudaContext& cu;
    const System& system;
    std::vector<int3> covalentFlagValues;
    std::vector<int2> polarizationFlagValues;
    CudaArray* multipoleParticles;
    CudaArray* molecularDipoles;
    CudaArray* molecularQuadrupoles;
    CudaArray* labFrameDipoles;
    CudaArray* labFrameQuadrupoles;
    CudaArray* sphericalDipoles;
    CudaArray* sphericalQuadrupoles;
    CudaArray* fracDipoles;
    CudaArray* fracQuadrupoles;
    CudaArray* field;
    CudaArray* fieldPolar;
    CudaArray* inducedField;
    CudaArray* inducedFieldPolar;
    CudaArray* torque;
    CudaArray* dampingAndThole;
    CudaArray* inducedDipole;
    CudaArray* inducedDipolePolar;
    CudaArray* inducedDipoleErrors;
    CudaArray* prevDipoles;
    CudaArray* prevDipolesPolar;
    CudaArray* prevDipolesGk;
    CudaArray* prevDipolesGkPolar;
    CudaArray* prevErrors;
    CudaArray* diisMatrix;
    CudaArray* diisCoefficients;
    CudaArray* extrapolatedDipole;
    CudaArray* extrapolatedDipolePolar;
    CudaArray* extrapolatedDipoleGk;
    CudaArray* extrapolatedDipoleGkPolar;
    CudaArray* inducedDipoleFieldGradient;
    CudaArray* inducedDipoleFieldGradientPolar;
    CudaArray* inducedDipoleFieldGradientGk;
    CudaArray* inducedDipoleFieldGradientGkPolar;
    CudaArray* extrapolatedDipoleFieldGradient;
    CudaArray* extrapolatedDipoleFieldGradientPolar;
    CudaArray* extrapolatedDipoleFieldGradientGk;
    CudaArray* extrapolatedDipoleFieldGradientGkPolar;
    CudaArray* polarizability;
    CudaArray* covalentFlags;
    CudaArray* polarizationGroupFlags;
    CudaArray* pmeGrid;
    CudaArray* pmeBsplineModuliX;
    CudaArray* pmeBsplineModuliY;
    CudaArray* pmeBsplineModuliZ;
    CudaArray* pmeIgrid;
    CudaArray* pmePhi;
    CudaArray* pmePhid;
    CudaArray* pmePhip;
    CudaArray* pmePhidp;
    CudaArray* pmeCphi;
    CudaArray* pmeAtomRange;
    CudaArray* lastPositions;
    CudaSort* sort;
    cufftHandle fft;
    CUfunction computeMomentsKernel, recordInducedDipolesKernel, computeFixedFieldKernel, computeInducedFieldKernel, updateInducedFieldKernel, electrostaticsKernel, mapTorqueKernel;
    CUfunction pmeSpreadFixedMultipolesKernel, pmeSpreadInducedDipolesKernel, pmeFinishSpreadChargeKernel, pmeConvolutionKernel;
    CUfunction pmeFixedPotentialKernel, pmeInducedPotentialKernel, pmeFixedForceKernel, pmeInducedForceKernel, pmeRecordInducedFieldDipolesKernel, computePotentialKernel;
    CUfunction recordDIISDipolesKernel, buildMatrixKernel, solveMatrixKernel;
    CUfunction initExtrapolatedKernel, iterateExtrapolatedKernel, computeExtrapolatedKernel, addExtrapolatedGradientKernel;
    CUfunction pmeTransformMultipolesKernel, pmeTransformPotentialKernel;
    CUevent syncEvent;
    CudaCalcMPIDGeneralizedKirkwoodForceKernel* gkKernel;
    static const int PmeOrder = 5;
    static const int MaxPrevDIISDipoles = 20;
};

/**
 * This kernel is invoked by MPIDMultipoleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMPIDGeneralizedKirkwoodForceKernel : public CalcMPIDGeneralizedKirkwoodForceKernel {
public:
    CudaCalcMPIDGeneralizedKirkwoodForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMPIDGeneralizedKirkwoodForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MPIDMultipoleForce this kernel will be used for
     */
    void initialize(const System& system, const MPIDGeneralizedKirkwoodForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Perform the computation of Born radii.
     */
    void computeBornRadii();
    /**
     * Perform the final parts of the force/energy computation.
     */
    void finishComputation(CudaArray& torque, CudaArray& labFrameDipoles, CudaArray& labFrameQuadrupoles, CudaArray& inducedDipole, CudaArray& inducedDipolePolar, CudaArray& dampingAndThole, CudaArray& covalentFlags, CudaArray& polarizationGroupFlags);
    CudaArray* getBornRadii() {
        return bornRadii;
    }
    CudaArray* getField() {
        return field;
    }
    CudaArray* getInducedField() {
        return inducedField;
    }
    CudaArray* getInducedFieldPolar() {
        return inducedFieldPolar;
    }
    CudaArray* getInducedDipoles() {
        return inducedDipoleS;
    }
    CudaArray* getInducedDipolesPolar() {
        return inducedDipolePolarS;
    }
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MPIDGeneralizedKirkwoodForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MPIDGeneralizedKirkwoodForce& force);
private:
    class ForceInfo;
    CudaContext& cu;
    const System& system;
    bool includeSurfaceArea, hasInitializedKernels;
    int computeBornSumThreads, gkForceThreads, chainRuleThreads, ediffThreads;
    MPIDMultipoleForce::PolarizationType polarizationType;
    std::map<std::string, std::string> defines;
    CudaArray* params;
    CudaArray* bornSum;
    CudaArray* bornRadii;
    CudaArray* bornForce;
    CudaArray* field;
    CudaArray* inducedField;
    CudaArray* inducedFieldPolar;
    CudaArray* inducedDipoleS;
    CudaArray* inducedDipolePolarS;
    CUfunction computeBornSumKernel, reduceBornSumKernel, surfaceAreaKernel, gkForceKernel, chainRuleKernel, ediffKernel;
};

/**
 * This kernel is invoked to calculate the vdw forces acting on the system and the energy of the system.
 */
class CudaCalcMPIDVdwForceKernel : public CalcMPIDVdwForceKernel {
public:
    CudaCalcMPIDVdwForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMPIDVdwForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MPIDMultipoleForce this kernel will be used for
     */
    void initialize(const System& system, const MPIDVdwForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MPIDVdwForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MPIDVdwForce& force);
private:
    class ForceInfo;
    CudaContext& cu;
    const System& system;
    bool hasInitializedNonbonded;
    double dispersionCoefficient;
    CudaArray* sigmaEpsilon;
    CudaArray* bondReductionAtoms;
    CudaArray* bondReductionFactors;
    CudaArray* tempPosq;
    CudaArray* tempForces;
    CudaNonbondedUtilities* nonbonded;
    CUfunction prepareKernel, spreadKernel;
};

/**
 * This kernel is invoked to calculate the WCA dispersion forces acting on the system and the energy of the system.
 */
class CudaCalcMPIDWcaDispersionForceKernel : public CalcMPIDWcaDispersionForceKernel {
public:
    CudaCalcMPIDWcaDispersionForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMPIDWcaDispersionForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MPIDMultipoleForce this kernel will be used for
     */
    void initialize(const System& system, const MPIDWcaDispersionForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MPIDWcaDispersionForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MPIDWcaDispersionForce& force);
private:
    class ForceInfo;
    CudaContext& cu;
    const System& system;
    double totalMaximumDispersionEnergy;
    CudaArray* radiusEpsilon;
    CUfunction forceKernel;
};

} // namespace OpenMM

#endif /*MPID_OPENMM_CUDAKERNELS_H*/
