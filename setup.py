# -- import packages: ---------------------------------------------------------
import setuptools

# -- read requirements: -------------------------------------------------------
with open('requirements.txt') as f:
    requirements = f.read().strip()

# -- read version: ------------------------------------------------------------
with open('ABCParse/__version__.py') as v:
    exec(v.read())

# -- setup: -------------------------------------------------------------------
setuptools.setup(
    name="ABCParse",
    version=__version__,
    python_requires=">=3.7.0",
    author="Michael E. Vinyard",
    author_email="mvinyard.ai@gmail.com",
    url="https://github.com/mvinyard/ABCParse",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    description="A better base class to automatically parse local args.",
    packages=setuptools.find_packages(),
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
    ],
    license="MIT",
)
