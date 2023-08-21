from setuptools import setup

setup(name='xbTools',
      version='0.0.6',
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
      install_requires=['numpy','matplotlib','netCDF4','scipy','datetime','geopandas','shapely','fiona','pytest'],
      zip_safe=False)
