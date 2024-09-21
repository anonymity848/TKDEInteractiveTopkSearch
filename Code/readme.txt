Readme (Interactive Top-k Query)
=========================
This package contains all source codes for 
a. Algorithm 2DSEG 
	1. It works for the special case of ITQ
	2. The code is in folder 2D
b. Algorithm 2DPI+ 
	1. It works for the special case of ITQ. 
	2. The code is in folder 2D
c. Algorithm Unify 
	1. It works for the general case of ITQ.
	2. The code is in folder HD
d. Algorithm UH-Simplex 
	1. It is an adapted existing algorithm.
	2. The code is in folder UH.
e. Algorithm ActiveRanking 
	1. It is an adapted existing algorithm.
	2. The code is in folder ActiveRanking.
f. Algorithm Pref-Learning 
	1. It is an adapted existing algorithm.
	2. The code is in folder PreferenceLearning.
g. Algorithm UH-Random 
	1. It is an adapted existing algorithm.
	2. The code is in folder UH.

Make sure there is a folder called "input/", a folder called "output/" and a file called 
"config.txt" under the working directory.
They will be used for storing the input/output files, some intermediate results and the 
input parameters.

Usage Step
==========
a. Compilation
	mkdir build
	cd build
	cmake ..
	make

	You will need to install the GLPK package (for solving LPs) at first.
	See GLPK webpage <http://www.gnu.org/software/glpk/glpk.html>.
	Then update the path in CMakeLists.txt
		set(INC_DIR /usr/local/Cellar/glpk/5.0/include)
		set(LINK_DIR /usr/local/Cellar/glpk/5.0/lib)
	Update path "/usr/local/Cellar/glpk/5.0" to the path you install the GLPK package
	
b. Execution
	./run

c. Config
	The config file contains the input parameters (whose format will be described in Appendix A).

Example
=======
Sample input (input/4d.txt) are provided. The dataset is described by four attributes

Try: ./run



Appendix A. Format of Config File
------------------------------------
The format is: AlgorithmName DatasetName k s u[1] u[2] ... u[d]
AlgorithmName - the name of the algorithm
DatasetName - the name of the dataset
k - the parameter that determines the quality of the output
s - the parameter that decides the output size
u[1] - the first dimension of the user's utility vector
u[2] - the second dimension of the user's utility vector
...
u[d] - the d-th dimension of the user's utility vector
For example, you might see
-----------------------
Unify23  4d.txt  10  5  0.216857  0.293712  0.469882  0.0195493  
-----------------------

