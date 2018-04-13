#!/bin/bash

wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh -b -p $PWD/miniconda2
(
    PATH=$PWD/miniconda2/bin:$PATH
    conda update conda
    conda config --append channels conda-forge
    conda install conda-build cmake pkgconfig emacs
    conda create -y -n ozodev petsc clangdev
)
git clone https://github.com/syl20bnr/spacemacs ./.emacs.d
