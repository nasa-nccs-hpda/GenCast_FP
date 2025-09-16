#!/bin/bash
#
# Build python module.
#
export SETUPTOOLS_USE_DISTUTILS=1 

f2py -c -m eta2xprs_ eta2xprs.F90 only: eta2xprs : 
