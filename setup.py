from setuptools import setup
from setuptools import find_packages

setup(name='DeepFRI',
      version='1.1.0',
      description='Implementation of Deep Functional Residue Identification',
      author='Vladimir Gligorijevic',
      author_email='vgligorijevic@flatironinstitute.org',
      url='https://github.com/bioinf-mcb/DeepFRI',
      download_url='https://github.com/bioinf-mcb/DeepFRI',
      license='FlatironInstitute',
      package_data={'DeepFRI': ['README.md']},
      packages=find_packages())
