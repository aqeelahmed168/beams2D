from setuptools import setup,  find_packages
from distutils.core import setup

PYTHON_REQUIRES = ">=3.9"

INSTALL_REQUIRES = [
    'numpy>=1.22',
    'pandas>=1.3.5',
    'scipy>=1.8'
]

PACKAGES = find_packages()

# list all the sub-packges/dirs
"""PACKAGES = [
    'beams2D',
    'beams2D.src',
    'beams2D.src.geom',
    'beams2D.src.mesh',
    'beams2D.src.elements',
    'beams2D.src.assembleSystem',
    'beams2D.inputs',
]
"""

setup(name='beams2D',
      version='0.1',
      description='Basic structural dynamics using Beam Elements in 2D',
      author='Aqeel Ahmed',
      license='MIT',
      python_requires=PYTHON_REQUIRES,
      install_requires=INSTALL_REQUIRES,
      packages=PACKAGES,
      zip_safe=True)


