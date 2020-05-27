import os
from os.path import relpath, join
from setuptools import setup
import versioneer


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def find_package_data(data_root, package_root):
    files = []
    for root, dirnames, filenames in os.walk(data_root):
        for fn in filenames:
            files.append(relpath(join(root, fn), package_root))
    return files


setup(
    name="offPELE",
    author="Martí Municoy",
    author_email="martimunicoy@gmail.com",
    description=("Open Forcield Open Energy Project"),
    license="MIT",
    keywords="molecular mechanics, forcefield, potential energy",
    url="http://github.com/martimunicoy/offPELE",
    packages=['offPELE'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 1 - Planning",
        "Natural Language :: English",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix"
    ]
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    package_data={'offPELE': find_package_data(
        'offPELE')}
)
