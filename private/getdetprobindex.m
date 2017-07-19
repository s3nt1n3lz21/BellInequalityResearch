function index = getdetprobindex(ni,di,mi,d,m)
% function to calculate the correct index of the variable specified by (ni,di,mi) in the vector vars
    index = (ni-1)*d*m+(di-1)*m+mi;