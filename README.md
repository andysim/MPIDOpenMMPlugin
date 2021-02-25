# About

This plugin implements polarizable multipole electrostatics, primarily aimed at
supporting the [MPID](https://doi.org/10.1063/1.4984113) formulation of the
CHARMM Drude force field.  The code supports multipoles up to octopoles, as
well as induced dipoles that may be either isotropic or anisotropic.

# Installation

Conda is *strongly* recommended for managing the environment and dependencies;
after [downloading](https://conda.io/docs/download.html) anaconda (make sure
you use `bash` or `zsh`).

To install the 7.5 verion of OpenMM into its own Conda environment called `mpid`, run
``` bash
conda create -n mpid openmm=7.5 cudatoolkit=10.2 swig mdtraj -c conda-forge
```
Make sure you request the version of the CUDA toolkit supported on your
cluster.  This example uses GCC to build; the speed of the C++ compiler is
irrelevant, because the CUDA code is the only fast code available in this
plugin.  Although the reference platform will run, it is very slow and
designed for correctness.

Then to build the code:

``` bash
conda activate mpid

export OPENMM_INSTALL_DIR=~/anaconda3/envs/mpid

module load cuda/10.2
module load cmake
module load gcc/8.2

# From the MPIDOpenMMPlugin top level directory
mkdir build
cd build
CXX=g++ cmake .. -DCMAKE_INSTALL_PREFIX=$OPENMM_INSTALL_DIR -DPYTHON_EXECUTABLE=`which python` -DOPENMM_DIR=$OPENMM_INSTALL_DIR -DCMAKE_CXX_FLAGS='-std=c++11'
make -j 4
make test
make install
make PythonInstall
```

Before running the code, make sure you load the conda environment and all modules used for building.
