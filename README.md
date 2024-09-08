
ænet-PyTorch
=====================

If you make use of the ænet-PyTorch interface, please cite the following reference:

**J. Lopez-Zorrilla<super>*</super>**, X.M. Aretxabaleta, I.W. Yeu, I. Etxebarria, H. Manzano, **N. Artrith<super>*</super>**, ænet-PyTorch: A GPU-Supported Implementation for Machine Learning Atomic Potentials Training, J. Chem. Phys. **158**, 164105 (2023). DOI: https://doi.org/10.1063/5.0146803 **OpenAccess**

<super>*</super>Contact:  jon.lopezz@ehu.eus or n.artrith@uu.nl

## **ænet**

<span id="sec:about"></span>

The Atomic Energy NETwork (**ænet**) package (http://ann.atomistic.net) is a collection of tools for the construction and application of atomic interaction potentials based on artificial neural networks (ANN). ANN potentials generated with **ænet** can then be used in larger scale atomistic simulations.


## **ænet-PyTorch**

**ænet-PyTorch** is an extension of that code to allow GPU-support for the training process of **ænet**, substituting the `train.x` training step. It is enterily written in PyTorch and includes new features that the previous code did not: the ability to fit reference forces in addition to energies with GPU support. **ænet-PyTorch** is fully compatible with all the **ænet** tools: interfaces with LAMMPS and TINKER, and ASE.

**M.S. Chen<super>*</super>**, T. Morawietz, H. Mori, T.E. Markland, **N. Artrith<super>*</super>**, AENET-LAMMPS and AENET-TINKER: Interfaces for Accurate and Efficient Molecular Dynamics Simulations with Machine Learning Potentials, J. Chem. Phys. **155**, 074801 (2021). doi: https://doi.org/10.1063/5.0063880

## *Forked-version*
I forked the [aenet-Pytorch repository](https://github.com/atomisticnet/aenet-PyTorch). 

1. I modified generate.x code. Now 
`aenet_generate_MPI.x` generates descriptors with the use of the MPI. 

2. I added a subroutine `generate_subroutine_MPI(inFile,ionum)`. You can generate descriptors by calling `generate_subroutine_MPI`. Here, `infile` is an input file like `generate.in`. `ionum` is a device number (integer). 
3. I made CMakeLists to use cmake. 

This forked version is made by Yuki Nagai (Information technology center, the University of Tokyo). 


# Installation
You can use CMake to install the ænet. 
`$ cd aenet_modified/`
`$ mkdir build`
`$ cd build`
`$ cmake ..`

Then, you can install both L-BFGS-B library and the ænet. 

## From the original README
**The following instruction is copied from the [aenet-Pytorch repository](https://github.com/atomisticnet/aenet-PyTorch).**

# Installation

<span id="sec:installation"></span>

## Installation of **ænet**

The modified version of ænet can be installed the same way as ænet. See its documentation or any of the tutorials available in the ænet for a comprehensive guide on how to install it. In short, these are the steps to follow:

1.  Compile the L-BFGS-B library
      - Enter the directory “aenet_modified/lib”
        
        `$ cd aenet_modified/lib`
    - Adjust the compiler settings in the “Makefile” and compile the library with
        
        `$ make`
    
    The library file `liblbfgsb.a`, required for compiling **ænet**,  will be created.

2.  Compile the **ænet** package
    
      - Enter the directory “aenet_modified/src”
        
        `$ cd aenet_modified/src`
    
      - Compile the ænet source code with an approproiate `Makefile.XXX`
        
        `$ make -f makefiles/Makefile.XXX`
    
    The following executables will be generated in “./bin”:
    
      - `generate.x`: generate training sets from atomic structure files
      - `train.x`: train new neural network potentials
      - `predict.x`: use existing ANN potentials for energy/force prediction



## Installation of **ænet-PyTorch**

**ænet-PyTorch** is enterily written in Python, using the PyTorch framework. In the following the required packages will be listed. It is highly recommended to install all the packages in an isolated Python environment, using tools such as virtual-env (https://virtualenv.pypa.io/en/latest/) or a conda environment (https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) .

  - `Python`: 3.7 or higher
  - `Numpy`: 1.19 or higher
  - `PyTorch`: 1.10 or higher
  - `CUDA`: 10.2 or higher (optional for GPU support) 

We will assume that `CUDA` and Python 3.7 are already installed.

1.  Install PyTorch 1.10

      - Installation using pip with CUDA support

        `$ pip install torch==1.10.1+cu102 -f https://download.pytorch.org/whl/cu102/torch_stable.html`

        or for only CPU usage
    
        `$ pip install torch==1.10.1+cpu -f https://download.pytorch.org/whl/cpu/torch_stable.html`


2.  Install Numpy

      - Numpy should be automatically installed when installing PyTorch.


## Installation of the tools for **ænet-PyTorch**

Compile the tools needed to make **ænet-PyTorch** compatible with **ænet**

  - Enter the directory "pytorch-aenet/tools"

    `cd tools`

  - Compile the tools

    `make`

The following exacutables will be generated in the same "tools/" directory:

  - `trainbin2ASCII.x`: convert the output from generate.x (`XXX.train`) to a format readable by **ænet-PyTorch** (`XXX.train.ascii`).
  - `nnASCII2bin.x`: convert the **ænet-PyTorch** output files (`XXX.nn.ascii`) to the usual binary files (`XXX.nn`).


# Usage

<span id="sec:usage"></span>

**ænet-PyTorch** is used in multiple steps:

1. ænet's regular `generate.x` is used to featurize the atomic structures in the reference data set;
2. The tool `trainbin2ASCII.x` is used to convert the training set file generated by `generate.x` (`XXX.train`) to a format readable by **ænet-PyTorch** (`XXX.train.ascii`);
3. Training is performed using **ænet-PyTorch**;
4. The trained potentials (`XXX.nn.ascii`) are converted back to the binary format used by ænet with `nnASCII2bin.x` (`XXX.nn`);
5. The potential can then be used by ænet's `predict.x` and other compatible interfaces, such as **ænet-LAMMPS**

An explanation of the input parameters can be found in the documentation in the `doc/` directory. An example of the usage of the code can be found in the `example/` directory.
