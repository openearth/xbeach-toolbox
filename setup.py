from pathlib import Path
from setuptools import setup

setup(
    name="xbTools",
    version="1.0.2",
    description="Toolbox to analyse and setup xbeach models",
    url="https://github.com/openearth/xbeach-toolbox",
    author="Menno de Ridder",
    author_email="menno.deridder@deltares.nl",
    long_description=Path("README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    license="GNU General Public License v3.0",
    packages=["xbTools", "xbTools.grid", "xbTools.general"],
    package_data={"": ["*.json"]},
    python_requires=">=3.12, <4",
    install_requires=[
        "numpy",
        "matplotlib>=3.5.0",
        "netCDF4>=1.6.0",
        "scipy>=1.10.0",
        "geopandas>=0.13.0",
        "shapely>=2.0.0",
        "fiona>=1.9.4.post1",
    ],
    zip_safe=False,
)
