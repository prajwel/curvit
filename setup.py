import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="curvit", 
    version="0.1.3",
    author="Prajwel Joseph",
    author_email="prajwel.joseph@gmail.com",
    description="lightcurves from UVIT data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/prajwel/curvit",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.4',
    install_requires=['numpy',
                      'matplotlib',
                      'astropy',
                      'scipy'],
)
