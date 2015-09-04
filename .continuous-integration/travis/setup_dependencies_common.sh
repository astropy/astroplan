#!/bin/bash -x

# CONDA
conda create --yes -n test -c astropy-ci-extras python=$PYTHON_VERSION pip
source activate test

# EGG_INFO
if [[ $SETUP_CMD == egg_info ]]
then
  return  # no more dependencies needed
fi

# PEP8
if [[ $MAIN_CMD == pep8* ]]
then
  pip install pep8
  return  # no more dependencies needed
fi

# CORE DEPENDENCIES
conda install --yes pytest Cython jinja2 psutil pytz

# NUMPY
if [[ $NUMPY_VERSION == dev ]]
then
  pip install git+http://github.com/numpy/numpy.git
  export CONDA_INSTALL="conda install --yes python=$PYTHON_VERSION"
else
  conda install --yes numpy=$NUMPY_VERSION
  export CONDA_INSTALL="conda install --yes python=$PYTHON_VERSION numpy=$NUMPY_VERSION"
fi

# ASTROPY
if [[ $ASTROPY_VERSION == dev ]]
then
  pip install git+http://github.com/astropy/astropy.git
else
  $CONDA_INSTALL astropy=$ASTROPY_VERSION
fi

# Now set up shortcut to conda install command to make sure the Python and Numpy
# versions are always explicitly specified.

# OPTIONAL DEPENDENCIES
if $OPTIONAL_DEPS
then
  # Note: nose is required to run the matplotlib image comparison tests
  $CONDA_INSTALL matplotlib nose
  pip install pyephem pytest-mpl
fi

# DOCUMENTATION DEPENDENCIES
# build_sphinx needs sphinx as well as matplotlib and wcsaxes (for plot_directive). 
if [[ $SETUP_CMD == build_sphinx* ]]
then
  $CONDA_INSTALL Sphinx Pygments matplotlib
  pip install wcsaxes
fi

# COVERAGE DEPENDENCIES
if [[ $SETUP_CMD == 'test --remote-data -V --coverage' ]]
then
  pip install coverage coveralls;
fi
