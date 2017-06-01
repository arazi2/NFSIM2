This package is a modified version of the NFSIM package.
The modifications are made by Abolfazl Razi in collaboration with Dr. Posner lab.
Please note that this program is still under test and development and a more robust version will be published soon. 
Last update : 05-30-2017. 
For all questions please email me at : abolfazl.razi@nau.edu.

The original version of NFSIM is available via the following link: 

https://github.com/msneddon/nfsim

For more information about NFSIM and how to use it, you can find a user manual at:

https://github.com/msneddon/nfsim/blob/master/NFsim_manual_v1.11.pdf


The main modifications are as follows:

-The following files are added:
	a) NFUtil/setting.hh  
		This file include some macros in order to activate/inactivate debugging message and modifications.
		
		By setting RAZI_DEBUG to different values, different sets of debugging messages are shown.
		
		Please also note that you can activate/inactivate messages with switch -v in command line.
		
		By setting RHS_FUNC to 1: the right hand side (RHS) function feature is activated.
		
		By uncommenting FIX_BUG1: a part of the code to fix of a bug regarding local functions is activated.
		
		By uncommenting FIX_A: the propensity calculation is modified. See the note at the end of this file for more details.
			
	
The following is the list of the most important changes:

- The makefile is update to include the new files in compile process

- A bug is fixed regarding using local functions.

	Two exemplary files are included to demonstrate the bug.
		lfb1.bngl includes two lines labeled as bug1 and bug2 which break down when running with NFSIM.
		lfb1.bngl includes equivalent lines labeled as Nobug1 and Nobug2 which work properly.
	The modified version executes both files properly without any error.
	
- A new feature is added to the code that enables using right hand side local functions in the reaction rules.

	This is a critical addition and enables evaluating and testing functions on the reaction products. 
	With this functionality, we can let a reaction to fire only if certain conditions are satisfied. 
	The RHS function are typically composite functions that include parameters, fixed values and local functions
	The results of local functions are expected to be 0 or 1 which means that a certain condition is met or not, whereas the fixed and global parameters define the reaction rate.
	The program when firing a reaction with RHS function, first chooses reactants using the mappingset objects, then clones the selected reactants. 
	Then reaction then is applied to the cloned reactants and the conditions [defined by RHS function] are evaluated on the resulting reaction products.
	If the conditions are satisfied, then the reaction is applied to the previously chosen real reactants. In both cases, the cloned reactants and corresponding mapping objects are killed. 
	
	
Some minor comments and notes to help run the software:

-Install a c++ compiler, for instance MinGW, add the corresponding path to the system environment parameters [in wondows].

-activate -std=gnu++11 or -std=c++11 to fix compile errors

-simple_system1 is a model based on the built-in simple_test in the code [the bngl file enables running in a usual way]

-Use DebugStructas a backup for Debug as well as the required makefile components [to run without the need for IDE]

-lbf1 and lbf2 are two bngl files with a minor difference [One ok and one broken model when used by the the older versions of NFSIM]
[local_func_bug is the original file that includes both ok and problematic reactions]

-t1func1 and t1func2 are two models, where the first one include regular LHS function and the second one includes a reaction with RHS function

-Use src -cn -v -xml t1func2.xml to run the model. [for more details src -h or read the manual of teh original code in the link provided above]
 
-Note about Propensity calculation:

It seems that there is a problem in calculating propensities. In the original version, the propensity of a reaction is calculated as base rate * n1*n2*..nk, where ni is the number of the ith reactant.  
This is problematic and perhaps a better option would be to use baserate * min(n1,n2,....n_k) or baserate * (n1,n2,....n_k)^(1/k).
Consider two reactions: 

	X(p~0)->X(p~1) 0.1  
	
	X(y)+Y(x)->X(y!1).Y(x!1)   0.1
	
If the initial number of both X and Y species is 100, the propensities of two reactions with the current calculations are

	p1 = 0.2 * 100 = 50   and    p2=0.2 * 100 * 100=5000 > 100.
p2 is more than the available reactants of types X(y) and Y(x). 
It seems incorrect for the second reaction [one may expect a lower rate for the second reaction or a rate below (100 * 0.2=20)]
This issue needs further test and investigation, so I do not recommend to activate it by setting FIX_A in setting.hh and using the original version.
