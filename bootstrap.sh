#!/usr/bin/env bash

# Ubuntu 16.04 LTS

update-locale LANG=en_US.UTF-8

apt-get update
apt-get install -y ntpd synaptic git emacs gfortran makedepf90 gnuplot5 petsc-dev global exuberant-ctags python-pygments elpa-helm m4 python-sympy ncview cmake
