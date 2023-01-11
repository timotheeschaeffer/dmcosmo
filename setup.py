'''
Created on 23 November 2021
@author: Timothée Schaeffer
Setup script
'''

from setuptools import setup, find_packages
#from distutils.core import setup

setup(name='dmcosmo',
      version='0.0.1',
      author='Timothee Schaeffer',
      author_email='timothee.schaeffer@uzh.ch',
      package_dir = {'dmcosmo' : 'src'},
      packages=['dmcosmo'],
      package_data={'share':['*'],},
      install_requires=['numpy','scipy', 'joblib', 'tqdm']
)
