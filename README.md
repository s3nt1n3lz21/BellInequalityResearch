# BellInequalityResearch
Summer 2017 Research Project

Code to calculate the maximum classical value of any bell inequality of the class of scenarios (n,2,m). Note throughout the scenarios are denoted (n,d,m) in that order where n is the number of parties, d is the number of measurement outcomes and m is the number of measurement settings. This is instead of the usual notation of (n,m,d) since the deterministic probabilities are written D_n_i(d_i|m_i) and makes the code easier to read and understand.

TODOS:
Create main file to loop over calcdimandclassical bound for each state and calculate the other quantities
Make n,m and d constants
Only return 2 argument from loopdetprobs
Only return non NaN rows from calcclassicalbound
Write the Readme file
Add comments
Include validation of calcclassicalbound (size of coefflist)
Some method of including only one copy of the constants n,m and d and other parameters so there are not many copies in memory
Set up Git so i don't have to put in the password every time
Generalise to any number of measurements and outcomes for each party
Generalise the calculation of tensor products of projectors to a state of d qubits. (Create a function)
Returns s = 0 if corrcoefflist is all zeros
Create a method to clear the values up to the nuber of rows it is full, don't use the zeros function to reset the whole array of zeros when most of them already are zero.