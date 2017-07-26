# BellInequalityResearch
Summer 2017 Research Project

Code to calculate the maximum classical value of any bell inequality of the class of scenarios (n,d,m). Note throughout the scenarios are denoted (n,d,m) in that order where n is the number of parties, d is the number of measurement outcomes and m is the number of measurement settings. This is instead of the usual notation of (n,m,d) since the deterministic probabilities are written D_n_i(d_i|m_i) and makes the code easier to read and understand.

TODOS:
Some of the properties could have been made local variables to the function, but weren't because of the recursion and would have taken up a lot of memory. However, they can be turned into handle variables.

OPTIMIZATION:
- Make m,n and d constants.
- When resetting arrays, clear the rows up to the number of rows they are full, don't use the zeros function to reset the whole array of zeros when most of them are already zero.

FUNCTIONALITY:
- Generalise to any number of measurements and outcomes for each party.
    - Create a function to get the correct settings from the index. This will just be a series of for loops and the settings are the current values of the counter variables.
- Generalise the calculation of tensor products of projectors in main to a state of d qubits.
- Change the code so that you can ask it to only calculate smax (quicker, don't have to store the deterministic probabilities) or smax and dimension.
- In main, collate the good states with integer coefficients (or half integer etc).

BUGS/VALIDATION:
- Doesn't seem to work for n = 4 comparing to literature.
- Include a check of the correlator coefficient list to see if it is all zeros.
- Check the algorithm works for first index.
- Instead of checking whether a coefficient is zero check whether it is very close to zero e.g within 10^7, this may help deal with precision errors.
- Check magnitude of vector b
- Go through code with these states, using breakpoints.

OTHER:
- Tidy up main code.
- Write the readme explanation and usage.
- Look at handle variables.

PERSONAL:
- Set up Git so i don't have to put in the password every time.