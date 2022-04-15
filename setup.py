import os
from setuptools import setup

version_py = os.path.join(os.path.dirname(__file__), 'kasel', 'version.py')
version = open(version_py).read().strip().split(
    '=')[-1].replace('"', '').strip()
long_description = """
``kasel`` is a pipeline for working with TB sequence data
"""

HERE = os.path.dirname(__file__)

with open(os.path.join(HERE, "requirements.txt"), "r") as f:
    install_requires = [x.strip() for x in f.readlines()]

setup(
    name="kasel",
    author="Adam Witney",
    description='A toolset for working with TB sequence data',
    long_description=long_description,
	url="https://github.com/awitney/kasel",
    version=version,
    install_requires=install_requires,
    requires=['python (>=3.5)'],
    packages=['kasel'],
    package_data={'kasel': ['workflow/Snakefile','workflow/**/*']},
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'kasel=kasel.kasel_main:main',
        ],
    },
    author_email="awitney@sgul.ac.uk",
    classifiers=[
        'Development Status :: 1 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)
