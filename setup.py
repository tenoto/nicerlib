# -*- coding: utf-8 -*-

# Learn more: https://github.com/kennethreitz/setup.py
# Reference : https://qiita.com/Kensuke-Mitsuzawa/items/7717f823df5a30c27077

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='nicerlib',
    version='0.4.0',
    description='Python library for data analyses of the NICER X-ray observatory',
    long_description=readme,
    author='Teruaki Enoto',
    author_email='teruaki.enoto@gmail.com',
    install_requires=['numpy'],
    url='https://github.com/tenoto/nicerlib',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)
