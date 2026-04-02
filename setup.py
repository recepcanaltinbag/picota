#!/usr/bin/env python
"""
Setup script for PICOTA - Pipeline for Identification of Composite Transposons from Assembly graphs
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read long description from README
here = Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8") if (here / "README.md").exists() else ""

setup(
    name="picota",
    version="1.0.0-rc1",
    description="Pipeline for Identification of Composite Transposons from Assembly graphs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="PICOTA Contributors",
    author_email="",
    url="https://github.com/recepcanaltinbag/picota",
    license="MIT",
    keywords=["bioinformatics", "transposon", "assembly-graph", "composite-transposon", "antibiotic-resistance"],
    
    project_urls={
        "Bug Tracker": "https://github.com/recepcanaltinbag/picota/issues",
        "Documentation": "https://github.com/recepcanaltinbag/picota/blob/main/README.md",
        "Source Code": "https://github.com/recepcanaltinbag/picota",
    },
    
    packages=find_packages(include=["picota", "picota.*"]),
    python_requires=">=3.8",
    
    install_requires=[
        "biopython>=1.79",
        "pandas>=1.3",
        "pyyaml>=6.0",
        "tqdm>=4.0",
        "requests>=2.27",
    ],
    
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.12",
        ],
    },
    
    entry_points={
        "console_scripts": [
            "picota=picota.picota:main",
        ],
    },
    
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
