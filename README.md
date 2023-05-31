# fetch
Frechet-based Efficient Time-series Comparison Heuristic

To compute the Fréchet distance we rely on the open-source implementation of https://github.com/derohde/Fred.

One should start by downloading the latest version of our source code. If git is installed, this step can be implemented as follows: 

$ git clone https://github.com/MORE-EU/complex-event-detection --recursive 

C++14 standard should be used for compiling the code. On most Linux distributions, the GNU Compiler Collection (GCC) is commonly used. To install GCC with C++14 support, the user can open a terminal and run the following command:  

$sudo apt-get install g++ 

The user can verify that the installation was successful by running 

$g++ --version  

in the terminal and checking if the version displayed includes C++14 support. Our implementation uses OpenMP. OpenMP is a library used for parallel programming in shared-memory systems. It is typically bundled with C and C++ compilers. To verify that OpenMP is already installed, one can check the version of the compiler as above. If it is not installed, one can install it using appropriate package 

managers based on the operating system. For example, in Ubuntu, one can run the following command: 

$sudo apt-get install libomp-dev 

Our implementation also makes use of pybind11 library. Pybind11 is a lightweight header-only library that facilitates the creation of Python bindings for C++ code. It simplifies the process of exposing C++ functions and classes to Python. To install pybind11, one can use pip by running the following command: 

$pip install pybind11 

A Makefile is available and can be used to build and assemble the project. As it is already mentioned above the makefile produces two executables in ./bin/ The first executable is motifs which corresponds to /src/main.cpp. and the second is subann.so which corresponds to the library file that can be loaded and used in Python as an extension module.  Here are some basic commands used with make (in the directory containing the makefile): 

$make: This is the primary command used to invoke make. When run without any arguments, it looks for a Makefile in the current directory and executes the default target specified in the Makefile. 

$make clean: This target is often defined in a Makefile to remove any generated or intermediate files created during the build process. Running make clean cleans up the project by deleting these files, preparing it for a fresh build. 


To run the provided Jupyter notebook, Python 3.8.5+ is required. For a guide on how to install python, one can visit https://www.python.org/. In addition, jupyter-notebook must also be installed by executing the following command: 

$ pip install pip install notebook 

The provided Jupyter notebook is MotifsFrechet_toy_example.ipynb
