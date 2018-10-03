# About

This plugin implements polarizable multipole electrostatics, primarily aimed at
supporting the [MPID](https://doi.org/10.1063/1.4984113) formulation of the
CHARMM Drude force field.  The code supports multipoles up to octopoles, as
well as induced dipoles that may be either isotropic or anisotropic.

# Installation

Conda is *strongly* recommended for managing the environment and dependencies;
after [downloading](https://conda.io/docs/download.html) anaconda (make sure
you use `bash` or `zsh`).

To install the nightly build of OpenMM into its own Conda environment called `mpid`, run
``` bash
conda create -n mpid numpy openmm -c omnia/label/devcuda92 -c omnia
```

Then to build the code:

``` bash
conda activate mpid

export OPENMM_INSTALL_DIR=~/anaconda3/envs/mpid

module load cuda/9.2
module load cmake
module load intel/18.0.1

# From the MPIDOpenMMPlugin top level directory
mkdir build
cd build
CC=icc CXX=icpc CXXFLAGS="-std=c++11" cmake .. -DCMAKE_INSTALL_PREFIX=$OPENMM_INSTALL_DIR -DPYTHON_EXECUTABLE=`which python` -DOPENMM_DIR=$OPENMM_INSTALL_DIR
make -j 4 install
make test
make PythonInstall
```

Before running the code, make sure you load the conda environment and modules used for building:

``` bash
conda activate mpid
module load cuda/9.2
module load intel/18.0.1
```
