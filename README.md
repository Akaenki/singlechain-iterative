# singlechian-iterative
This code requires [OpenBLAS](https://github.com/xianyi/OpenBLAS). [OpenMP](https://computing.llnl.gov/tutorials/openMP/) is optional but recommanded.
This code can run iterative method of dilute system
I only test it with gcc, may not be compatible with icc or clang. The code is *only* compativble with C99. 

# Usage
I use Makefile to compile the code, Some common input is also in the Makefile.
##Makefile
* `CC=gcc-5` change this to the compiler you are using I am using gcc/5.3 here.
* `OMPGLAG=-fopenmp` the OpenMP flag, for icc user use `-qopenmp`. Clang does not support OpenMP therefore leave it blank
* `HPATH=` and `LIBPARH=` theses are the library and header path for OpenBLAS. See it's [wiki](https://computing.llnl.gov/tutorials/openMP/) for detail.

* Other Inputs: there are some inputs in the Makefile, the meaning of them are described in it. 
* `TAEGETS=` the target program name, make sure to change the name before its rule when changing the target name.

To compile the code just run the following in terminal:
````
make clean && make
````

The inputs that are not included in the Makefile need to be specified through arguments
## Input arguments
* Running mode use flog `-m`, 0 - FD, 1 - HI. Default 0. 
* Chain length use flag `-n`. Default value 20
* Flow type use flag `-f`, 0 - no flow, 1 - planar elongational, 2 - steady shear. Default value 0
* Flow rate use flag `-r`, **must be specified** if flow type != 0. 
* Number of OpenMP threads use flag `-t`. Default 1 (no parallelism)

Example:
````
./program -t 4 -n 100 -f 1 -r 0.01 -m 1
````
# Others
* To automatically start the next iteration after the first iteration finished I included one method using the `run.sh` file. 

* The code can also calculate (1) radius of gyration (2) end-to-end vectors (3) self-diffusivity (4) extension. To use them uncomment the functions in `main.c`

* The code can also simulate ring polymers, just replace the right functions. For more detail see the comments in code.  
