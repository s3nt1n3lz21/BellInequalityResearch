Title: Bell Inequality Classical Bound And Dimension Calculator.
Author: Neil Smith
Date: 26/07/2017
Contact Info: neilcsmith21@gmail.com
Affiliation: The University Of Sheffield

DESCRIPTION:
An algorithm to calculate the maximum classical value and dimension of any bell inequality of the scenarios (n,d,m) where n is the number of parties, d is the number of measurement outcomes and m is the number of measurement settings, here it is assumed that these numbers are the same for each party.
To ensure it always gets the correct answer, it simply iterates over all possible combinations of the single party deterministic probabilities (0 or 1) which are written as D_ni(di,mi) for party ni with measurement outcome di and setting mi. This can be thought of as changing a binary number from all zeros 000... to all ones 111... Since there are nmd of these deterministic probabilities or variables then there will be 2^(nmd) possible states to explore, the problem is therefore exponential in time but can be made polynomial in memory. 
The simplest method is to just write out a series of for loops, but in general there will be nmd for loops and so to generalise the process, recursion is used. A function loops over the possible values of a variable then calls the function again with one less variable to loop over, if there are no more then it performs the calculation. However doing this with a normal function would lead to large numbers of copies of variables in memory, instead a class design is used where the data the algorithm has to perform calculations upon is made a property so any part of the inner nested recursive functions can access the data at any time and ensures there is only one copy of the data.

USAGE:
To use the algorithm simply specify the scenario and the correlator coefficient list associated with the Bell Inequality when calling the calcdimandclassicalbound function. The list should be a row vector (m+1)^n long. 
As an example of the form of the input consider the CHSH inequality in the scenario (2,2,2), the list has the form [0 0 0 0 1 1 0 1 -1]. The first element represents all the parties making no measurements, a party making no measurement represents a measurement setting of 0 and so this would be written as 00. The second represents the second person performing the first measurement so 01 and the third represents the second person performing their second measurement so 02. The next would be 10 and so on, representing the incrementation of a base n number.
So this represents a Bell Inequality of the form 0*[00]+0*[01]+0*[02]+0*[10]+1*[11]+1*[12]+0*[20]+1*[21]-1*[22].

FILELIST:
- "cbanddimcalc.m" The algorithm implemented as a class.
- "calcdimandclassicalbound.m" The function that acts as a wrapper for the class, this is what a user would call.
- "main.m" An 

#######################################################

TODOS:

OPTIMIZATION:
- Make m,n and d constants.
- When resetting arrays, clear the rows up to the number of rows they are full, don't use the zeros function to reset the whole array of zeros when most of them are already zero.
- Maybe have a better implementation using handle variables, is this possible?

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
- Instead of checking whether a coefficient is zero check whether it is very close to zero e.g within 10^7, this may help deal with precision errors and speed up the calculation.