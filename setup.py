#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""atomic_charges_step
A SEAMM plug-in for computing atomic (partial) charges from a converged
electron density, using real-space density-partitioning programs such as
DDEC6 (Chargemol) and Bader (Henkelman).
"""
import sys
from setuptools import setup, find_packages
import versioneer


# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

with open('requirements.txt') as fd:
    requirements = fd.read()

setup(
    name='atomic_charges_step',
    author="Paul Saxe",
    author_email='psaxe@molssi.org',
    description=__doc__.splitlines()[1],
    long_description=readme + '\n\n' + history,
    long_description_content_type='text/x-rst',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    url='https://github.com/molssi-seamm/atomic_charges_step',

    packages=find_packages(include=['atomic_charges_step']),

    include_package_data=True,

    setup_requires=[] + pytest_runner,

    install_requires=requirements,

    test_suite='tests',

    platforms=['Linux',
               'Mac OS-X',
               'Unix',
               'Windows'],

    zip_safe=True,

    keywords=['SEAMM', 'SEAMMplugin', 'flowchart'],
    classifiers=[
        'Environment :: Plugins',
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    entry_points={
        'console_scripts': [
            'atomic-charges-step-installer='
            'atomic_charges_step.__main__:run',
        ],
        'org.molssi.seamm': [
            'Atomic Charges = atomic_charges_step:AtomicChargesStep',
        ],
        'org.molssi.seamm.tk': [
            'Atomic Charges = atomic_charges_step:AtomicChargesStep',
        ],
    }
)
