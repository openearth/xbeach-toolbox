from setuptools import setup

setup(name='xbTools',
      version='1.0.0',
      description='Toolbox to analyse and setup xbeach models',
      url='https://github.com/openearth/xbeach-toolbox',
      # Authors: Menno de Ridder (menno.deridder@deltares.nl), Cas van Bemmelen (cas.van.bemmelen@witteveenbos.com)
      # PyPI does not have an easy way to specify multiple authors.
      author='Menno de Ridder',
      author_email='menno.deridder@deltares.nl',
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown',
      license='GNU General Public License v3.0',
      packages=['xbTools','xbTools.grid','xbTools.general'],
      package_data={'': ['*.json']}, # always include .json files in the package contents
      python_requires='>=3, <4',
      install_requires=['numpy','matplotlib>=3.5.0','netCDF4>=1.6.0','scipy>=1.10.0','datetime','geopandas>=0.13.0','shapely>=2.0.0','fiona>=1.9.4.post1','pytest'],
      zip_safe=False)
