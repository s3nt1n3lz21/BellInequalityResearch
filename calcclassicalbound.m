% Calculate the classical bound of the expression by calculating the expression for all possible values of the single party deterministic probabilities (vars)
function [smax,detprobsgivesmax,numdetprobs] = calcclassicalbound(n,d,m,corrcoefflist)
detprobs = zeros(1,n*d*m);
[finalsmax,finaldetprobsgivesmax,finalamountdetprobsgivesmax] = loopdetprobs(n*d*m,detprobs,corrcoefflist,n,d,m);
smax = finalsmax;
detprobsgivesmax = finaldetprobsgivesmax;
numdetprobs = finalamountdetprobsgivesmax;