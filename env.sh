#!/bin/bash

thisdir=$(dirname $(readlink -f ${BASH_SOURCE[0]}))

export PATH=$thisdir/miniconda3/bin:$PATH
export PKG_CONFIG_PATH=$thisdir/miniconda3/lib/pkgconfig:$PKG_CONFIG_PATH
alias emacs="HOME=$thisdir emacs"
