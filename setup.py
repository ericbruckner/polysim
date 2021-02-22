from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="poly_dpd", 
    version="0.1.0",
    author="Eric Bruckner",
    author_email="eric.p.bruckner@gmail.com",
    description="A library for running DPD simulations of simple polymers with HOOMD-blue",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(include = 'poly_dpd'),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent"
    ],
    python_requires='>=3',
)
