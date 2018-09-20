# About

This plugin implements polarizable multipole electrostatics, primarily aimed at
supporting the [MPID](https://doi.org/10.1063/1.4984113) formulation of the
CHARMM Drude force field.  The code supports multipoles up to octopoles, as
well as induced dipoles that may be either isotropic or anisotropic.

# Installation

This plugin requires a specially modified version of OpenMM, which currently
lives on the [mpid
branch](https://git.lobos.nih.gov/andysim/OpenMM/commits/mpid) in the LCB's
GitLab repo.

##### Prerequisites
Conda is *strongly* recommended for managing the environment and dependencies;
after [downloading](https://conda.io/docs/download.html) anaconda (make sure
you use `bash` or `zsh`), you should be able to create an environment:

```
conda create -n mpid python=3.6 numpy
```
##### OpenMM
To install OpenMM, make sure you check out the `mpid` branch of the above repo
and from a build subdirectory, call CMake before building and installing.  For
example, on the NIH LoBoS cluster, the following should be run:
``` bash
# Specify where to install the mpid version of OpenMM 
export OPENMM_INSTALL_DIR=/path/to/locatation/where/openmm/and/mpid/plugin/should/live

module load cuda/9.2
module load cmake
module load intel/18.0.1
module load doxygen
module load swig/3.0.12

# Load the mpid conda environment, which sets the path up correctly
source activate mpid
# Call CMake to make sure that it's all
CC=icc CXX=icpc CXXFLAGS="-std=c++11" cmake .. -DCMAKE_INSTALL_PREFIX=$OPENMM_INSTALL_DIR -DPYTHON_EXECUTABLE=`which python`

# Run the build (set the number to how every many cores you have)
make -j 24
make test
make install
make PythonInstall
```

##### The plugin
``` bash
CC=icc CXX=icpc CXXFLAGS="-std=c++11" cmake .. -DCMAKE_INSTALL_PREFIX=$OPENMM_INSTALL_DIR -DPYTHON_EXECUTABLE=`which python` -DOPENMM_DIR=$OPENMM_INSTALL_DIR
make -j 24
make test
make PythonInstall
```

##### Running from python
``` bash
# Before running, the library locations need to be setup (don't forget to load the approprate Intel and CUDA modules, also)
export OPENMM_LIB_PATH=${OPENMM_INSTALL_DIR}/lib
export OPENMM_PLUGIN_DIR=${OPENMM_LIB_PATH}/plugins
export LD_LIBRARY_PATH=${OPENMM_LIB_PATH}:$LD_LIBRARY_PATH #N.B. this is $DYLD_LIBRARY_PATH on macOS instead!
```
