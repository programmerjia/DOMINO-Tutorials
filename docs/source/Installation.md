# Installation
The domino-spatial package is developed based on the pytorch framework and can be implemented on both GPU and CPU. We recommend running the package on GPU.  

For convenience, we suggest using a separate conda environment for running DOMINO. Please ensure anaconda3 and the required environment is installed.

Create conda environment and install domino-spatial package.

```bash
# create an environment called DOMINO

conda create -n DOMINO python=3.8

# activate your environment

conda activate DOMINO

# install the required environment

# install package

pip install domino-spatial
```

Then test whether the domino-spatial package has been installed successfully:

```bash
# test the installation

domino-spatial --version
```

If the package name and version number are output, the installation is successful.
