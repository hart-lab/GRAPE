from setuptools import setup, find_packages
 
from grape.utils.version import __version__
 

setup(name='GRAPE', version=__version__,
      author='Juihsuan Chou, Traver Hart',
      author_email='',
      description='GRAPE: Genetic interaction Regression Analysis of Pairwise Effects',
      license='MIT', packages=find_packages(),
      entry_points={'console_scripts': ['grape = grape.cli:__main__']},
      install_requires=['numpy==1.26.2', 'pandas==2.1.4', 'scipy==1.11.4', 'statsmodels==0.14.1',
                        'scikit-learn==1.3.2'])
