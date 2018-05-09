# nicerlib (user script library for NICER)
Python library for user analyses pipeline and scripts of the NICER X-ray observatory.

The current version requires NICERDAS version 2018-03-01_V003 included HEASoft version 6.24.

## Description

To use these scripts add *basedir*/nicerlib/scripts to your PATH and add *basedir*/nicerlib to your PYTHONPATH, where *basedir* is wherever you cloned nicerlib. After modifying environment setups of *basedir*/setenv/nicerlib_setenv.bashrc, you may run...

```
source basedir/setenv/nicerlib_setenv.bashrc
```


## Preparation 
You can skip this step if you have already installed required python library. It is recommended to construct nicerlib setups under a'pyenv' environment to avoid potential conflicts among python library versions. Followings are an installation example of 'pyenv'.

```
git clone git://github.com/yyuu/pyenv.git ~/.pyenv

echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.zshrc
echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.zshrc
echo 'eval "$(pyenv init -)"' >> ~/.zshrc

pyenv install 2.7.14
pyenv rehash
 
pyenv global 2.7.14
pyenv versions 
which python
pyenv which python2.7
```

### (optional) pyenv virtualenv 
```
git clone https://github.com/yyuu/pyenv-virtualenv.git ~/.pyenv/plugins/pyenv-virtualenv
```

add the following line in the .zshrc (or equivalent).
```
eval "$(pyenv virtualenv-init -)"
```

### Required software and python library 
- [HEASoft](https://heasarc.nasa.gov/lheasoft/)
- astropy (for time conversion only)
Required python modules can be installed for example...

```
pip install --upgrade pip

pip install astropy
pip show astropy
```

### (optional) PINT
PINT[[http://nanograv-pint.readthedocs.io/en/latest/index.html]]
see the installation [[http://nanograv-pint.readthedocs.io/en/latest/installation.html#prerequisites]]


## Installation
You can download the code as follows. 
```
git clone https://github.com/tenoto/nicerlib.git
```

## How to use
TBD 

### Directory

- nicerlib : a set of libraries for NICER data analyses. 
- scripts : command-line scripts. some of them are based on HEASoft.
- setenv : code for environment setups.
- test : test scripts and example codes to run the library in command line. 

### Scripts


## Change History 
* 2018-04-02 : niprefilter2 is included. 
* 2018-03-19 : First nicerlib script 

## Reference:
- (in Japanese) https://qiita.com/Kensuke-Mitsuzawa/items/7717f823df5a30c27077 

## Quick view of the user standard process of NICER data 
- nicercal (calibration)
- nimaketime (time screening)
- nicermergeclean (merge and clean event data)
- or... nicerl2 (calibration + screening) = nicercal + nimaketime + nicermergeclean 
- xselect (for extracting spectra an light curves)
- barycorr (barycentric correction)





