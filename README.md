# About

This plugin implements polarizable multipole electrostatics, primarily aimed at supporting the [MPID](https://doi.org/10.1063/1.4984113) formulation of the CHARMM Drude force field.  The code supports multipoles up to octopoles, as well as induced dipoles that may be either isotropic or anisotropic.

# Installation

This plugin requires a specially modified version of OpenMM, which currently lives on the [mpid branch](https://git.lobos.nih.gov/andysim/OpenMM/commits/mpid) in the LCB's GitLab repo.
##### Prerequisites
Conda is *strongly* recommended for managing the environment and dependencies; after [downloading](https://conda.io/docs/download.html) anaconda (make sure you use `bash` or `zsh`), you should be able to create an environment:

```
conda create -n mpid python=3.6 numpy
```
##### OpenMM
To install OpenMM, make sure you check out the `mpid` branch of the above repo and from a build subdirectory, run
``` bash
# Specify where to install the mpid version of OpenMM 
export OPENMM_INSTALL_DIR=/path/to/locatation/where/openmm/should/live
# Load the mpid conda environment, which sets the path up correctly
source activate mpid
# Call CMake to make sure that it's all
CXXFLAGS="-std=c++11" cmake .. -DCMAKE_INSTALL_PREFIX=$OPENMM_INSTALL_DIR -DPYTHON_EXECUTABLE=`which python`
# Run the build (set the number to how every many cores you have)
make -j 8 
make PythonInstall
```

##### The plugin
# Specify where to install the MPID plugin
``` bash
export PLUGIN_INSTALL_DIR=/path/to/locatation/where/the/mpid/plugin/should/live
CXXFLAGS="-std=c++11" cmake .. -DOPENMM_DIR=$OPENMM_INSTALL_DIR -DCMAKE_INSTALL_PREFIX=$PLUGIN_INSTALL_DIR
make -j 8
make PythonInstall
```
##### Running from python
export OPENMM_PLUGIN_DIR=$PLUGIN_INSTALL_DIR
export $DYLD_LIBRARY_PATH=${PLUGIN_INSTALL_DIR}/lib:$DYLD_LIBRARY_PATH #macOS
export $LD_LIBRARY_PATH=${PLUGIN_INSTALL_DIR}/lib:$LD_LIBRARY_PATH     #linux

