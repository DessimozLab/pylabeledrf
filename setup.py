import setuptools
import sys

name = 'pylabeledrf'
__version__ = "Undefined"
for line in open('{}/__init__.py'.format(name)):
    if (line.startswith('__version__')):
        exec(line.strip())

print(setuptools.find_packages())


with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
     name=name,  
     version=__version__,
     author="Christophe Dessimoz",
     author_email="christophe.dessimoz@unil.ch",
     description="A package to compute the Robinson Fould distance extended to labeled topologies.",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/dessimozlab/pylabeledrf",
     download_url="https://github.com/DessimozLab/pylabeledrf/archive/v0.2.2.tar.gz",
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
     install_requires=[
          'dendropy>=4',
     ],
 )
