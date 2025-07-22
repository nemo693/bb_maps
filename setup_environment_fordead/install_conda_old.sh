#!/bin/bash
set -x
CONDA_DIR="/opt/conda"
MINICONDA_VERSION="latest"
cp eurac_env.yml /opt/
cd /opt/
mkdir -p $CONDA_DIR &&                                                                              \
wget --quiet https://repo.continuum.io/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh && \
/bin/bash Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -f -b -p $CONDA_DIR &&                    \
rm Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh &&                                               \
$CONDA_DIR/bin/conda config --prepend channels conda-forge --system &&                               \
$CONDA_DIR/bin/conda config --prepend channels conda-forge/label/dev &&                              \
$CONDA_DIR/bin/conda update --all -y &&                                                              \
$CONDA_DIR/bin/conda env create  --file eurac_env.yml &&                           \
$CONDA_DIR/bin/conda clean -tipsy && \
$CONDA_DIR/bin/conda config --set auto_activate_base false
$CONDA_DIR/bin/conda init

