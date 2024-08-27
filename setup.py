from setuptools import setup, find_packages

setup(
    name="miasort",
    version="0.1.2",
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
)