# nicerlib
Python library for analyses pipeline and scripts of the NICER X-ray observatory.

## Description
To use these scripts add <basedir>/nicerlib/scripts to your PATH and add <basedir>/nicerlib to your PYTHONPATH, where <basedir> is wherever you cloned nicerlib.

## Preparation 
It is recommended to construct nicerlib setups under a'pyenv' environment to avoid potential conflicts among python library versions. Followings are an installation example of 'pyenv'.

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

Required python modules can be installed for example...

```
pip install --upgrade pip

pip install astropy
pip show astropy
```

## Installation
You can download the code as follows. 
```
git clone https://github.com/tenoto/nicerlib.git
```

## How to use
TBD 


## Change History 
* 2018-04-02 : niprefilter2 is included. 
* 2018-03-19 : First nicerlib script 

## Reference:
- (in Japanese) https://qiita.com/Kensuke-Mitsuzawa/items/7717f823df5a30c27077 
