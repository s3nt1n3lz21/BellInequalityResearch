# BellInequalityResearch
Summer 2017 Research Project

Code to calculate the maximum classical value of any bell inequality of the class of scenarios (n,2,m). Note throughout the scenarios are denoted (n,d,m) in that order where n is the number of parties, d is the number of measurement outcomes and m is the number of measurement settings. This is instead of the usual notation of (n,m,d) since the deterministic probabilities are written D_n_i(d_i|m_i) and makes the code easier to read and understand.

TODOS:
Check second 4,2,2 expression
Make n,m and d constants
Write the Readme file
Add comments
Include validation of calcclassicalbound (size of coefflist)
Set up Git so i don't have to put in the password every time
Generalise to any number of measurements and outcomes for each party
Generalise the calculation of tensor products of projectors to a state of d qubits.
Returns s = 0 if corrcoefflist is all zeros, Just skip the first index
Create a method to clear the values up to the nuber of rows it is full, don't use the zeros function to reset the whole array of zeros when most of them already are zero.
Create a getsettings function that is just a series of for loops and the settings are the current values of the counter variables, generalise it.
Change Calcprobdist to a class
Change the code so that you can ask it to only calculate smax (quicker) or smax and dimension.

Notes
Doesn't seem to work for n = 4 comparing to literature