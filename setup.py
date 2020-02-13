import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="micropy",
    version="0.6",
    author="Matteo Pedrotti",
    author_email="matteo.pedrotti@strath.ac.uk",
    description="Tool to investigate pore space",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/PedroMat8/pore_sizer",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
