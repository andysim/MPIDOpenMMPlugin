/* -------------------------------------------------------------------------- *
 *                              OpenMMMPID                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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

#include "MPIDCudaKernelFactory.h"
#include "MPIDCudaKernels.h"
#include "CudaPlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/windowsExport.h"

using namespace OpenMM;

#ifdef OPENMM_BUILDING_STATIC_LIBRARY
static void registerPlatforms() {
#else
extern "C" OPENMM_EXPORT void registerPlatforms() {
#endif
}

#ifdef OPENMM_BUILDING_STATIC_LIBRARY
static void registerKernelFactories() {
#else
extern "C" OPENMM_EXPORT void registerKernelFactories() {
#endif
    try {
        Platform& platform = Platform::getPlatformByName("CUDA");
        MPIDCudaKernelFactory* factory = new MPIDCudaKernelFactory();
        platform.registerKernelFactory(CalcMPIDBondForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcMPIDAngleForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcMPIDInPlaneAngleForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcMPIDPiTorsionForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcMPIDStretchBendForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcMPIDOutOfPlaneBendForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcMPIDTorsionTorsionForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcMPIDMultipoleForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcMPIDGeneralizedKirkwoodForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcMPIDVdwForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcMPIDWcaDispersionForceKernel::Name(), factory);
    }
    catch (...) {
        // Ignore.  The CUDA platform isn't available.
    }
}

extern "C" OPENMM_EXPORT void registerMPIDCudaKernelFactories() {
    try {
        Platform::getPlatformByName("CUDA");
    }
    catch (...) {
        Platform::registerPlatform(new CudaPlatform());
    }
    registerKernelFactories();
}

KernelImpl* MPIDCudaKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    CudaPlatform::PlatformData& data = *static_cast<CudaPlatform::PlatformData*>(context.getPlatformData());
    CudaContext& cu = *data.contexts[0];

    if (name == CalcMPIDBondForceKernel::Name())
        return new CudaCalcMPIDBondForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcMPIDAngleForceKernel::Name())
        return new CudaCalcMPIDAngleForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcMPIDInPlaneAngleForceKernel::Name())
        return new CudaCalcMPIDInPlaneAngleForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcMPIDPiTorsionForceKernel::Name())
        return new CudaCalcMPIDPiTorsionForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcMPIDStretchBendForceKernel::Name())
        return new CudaCalcMPIDStretchBendForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcMPIDOutOfPlaneBendForceKernel::Name())
        return new CudaCalcMPIDOutOfPlaneBendForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcMPIDTorsionTorsionForceKernel::Name())
        return new CudaCalcMPIDTorsionTorsionForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcMPIDMultipoleForceKernel::Name())
        return new CudaCalcMPIDMultipoleForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcMPIDGeneralizedKirkwoodForceKernel::Name())
        return new CudaCalcMPIDGeneralizedKirkwoodForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcMPIDVdwForceKernel::Name())
        return new CudaCalcMPIDVdwForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcMPIDWcaDispersionForceKernel::Name())
        return new CudaCalcMPIDWcaDispersionForceKernel(name, platform, cu, context.getSystem());

    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
