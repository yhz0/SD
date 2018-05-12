# sd

Prerequisite: CPLEX should be available on your machine. (If not, follow this link to get a trial version: http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/).

1). Download and compile SD with the following command. Simply copy and paste the following to your terminal (only support unix-like systems such as OS X and Ubuntu):

`git clone https://github.com/imliuyifan/sd && cd sd/src && sudo make && cd ../instance && ln -s ../src/sd`

(note: In order to locate the path to cplex header and library file, "sudo make" is required as you can see. You can check the content of the makefile to verify this.)

2). Then excute sd by typing the following into your terminal:

`./sd`

At the prompt, enter an instance name for example: `pgp2`

(note: all instances are stored in the sdinput folder.)

Results will be stored in sdoutput/pgp2. Please check pgp2.detailed_soln.out for solutions and time_sample.out for CPU time and number of samples used.
