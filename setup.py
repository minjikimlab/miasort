from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension(
        "miasort.util",
        ["miasort/util.pyx"],
        include_dirs=[numpy.get_include()],  # Include NumPy headers
    ),
    Extension(
        "miasort.process_complex_cython",  # Name of the compiled module
        ["miasort/process_complex_cython.pyx"],  # Path to the .pyx file
        include_dirs=[numpy.get_include()],  # Include NumPy headers if needed
    )
]

setup(
    name="miasort",
    version="0.1.6",
    author="Zichen Zhang",
    author_email="zhangzzc@umich.edu",
    description="A Tool for Multiplex Chromatin Interaction Analysis by Efficiently Sorting Chromatin Complexes",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/minjikimlab/mia-sort",
    packages=find_packages(),
    install_requires=[
        'matplotlib==3.9.2',
        'pybedtools==0.10.0',
        'setuptools==70.3.0',
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.11',
    ext_modules=cythonize(extensions)
)

