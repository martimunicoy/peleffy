import os
from os.path import relpath, join
from setuptools import setup
import versioneer


def find_package_data(data_root, package_root):
    files = []
    for root, dirnames, filenames in os.walk(data_root):
        for fn in filenames:
            files.append(relpath(join(root, fn), package_root))
    return files


setup(
    name="peleffy",
    author="Martí Municoy",
    author_email="martimunicoy@gmail.com",
    description=("PELE Force Field Yielder"),
    license="MIT",
    keywords="molecular mechanics, forcefield, potential energy",
    url="http://github.com/martimunicoy/peleffy",
    packages=[
        'peleffy',
        'peleffy/tests',
        'peleffy/data',
        'peleffy/template',
        'peleffy/solvent',
        'peleffy/topology',
        'peleffy/forcefield',
        'peleffy/utils',
        'peleffy/charge',
    ],
    long_description='The offpele (Open Force Field to PELE) is ' \
        + 'a Python package that builds PELE-compatible force ' \
        + 'field templates.',
    classifiers=[
        "Development Status :: 1 - Planning",
        "Natural Language :: English",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix"
    ],
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    package_data={'peleffy': find_package_data(
        'peleffy/data', 'peleffy')},
)
