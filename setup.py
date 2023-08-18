from setuptools import setup

setup(name='xbeach-tools-bemc',
      version='0.0.2',
      description='Toolbox to analyse and setup xbeach models',
      url='https://github.com/openearth/xbeach-toolbox',
      author='Menno de ridder',
      author_email='menno.deridder@deltares.nl',
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown',
      license='MIT',
      packages=['xbTools'],
      python_requires='>=3, <4',
      install_requires=['numpy','matplotlib','netCDF4','scipy','datetime','geopandas','shapely','fiona','pytest'],
      zip_safe=False)
