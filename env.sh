#!/bin/bash

thisdir=$(dirname $(readlink -f ${BASH_SOURCE[0]}))

source $thisdir/miniconda2/etc/profile.d/conda.sh
conda activate ozodev
#export PKG_CONFIG_PATH=$thisdir/miniconda3/lib/pkgconfig:$PKG_CONFIG_PATH
alias emacs="HOME=$thisdir emacs"
