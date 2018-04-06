# nicerlib
Python script for data analyses of the NICER X-ray observatory.

reference:
- (in Japanese) https://qiita.com/Kensuke-Mitsuzawa/items/7717f823df5a30c27077 

To use these scripts add <basedir>/nicerlib/scripts to your PATH and add <basedir>/nicerlib to your PYTHONPATH, where <basedir> is wherever you cloned nicerlib.


git clone git://github.com/yyuu/pyenv.git ~/.pyenv

echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.zshrc
echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.zshrc
echo 'eval "$(pyenv init -)"' >> ~/.zshrc

pyenv install 2.7.14
pyenv rehash

pyenv global 2.7.14
pyenv versions 
which python


pip install --upgrade pip

pip install astropy
pip show astropy


git clone https://github.com/tenoto/nicerlib.git


# Change History 
* 2018-04-02 : niprefilter2 is included. 
* 2018-03-19 : First nicerlib script 