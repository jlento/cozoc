#!/bin/bash

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -p $PWD/miniconda3
(
    PATH=$PWD/miniconda3/bin:$PATH
    conda install -c spectraldns netcdf4-parallel
    conda install -c conda-forge petsc emacs
)
git clone https://github.com/syl20bnr/spacemacs ./.emacs.d
