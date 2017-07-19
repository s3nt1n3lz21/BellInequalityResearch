function [dimension,smax] = calcdimandclassicalbound(n,d,m,corrcoefflist)
% CALCDIMANDCLASSICALBOUND Calculates the classical bound of a bell
% inequality and it's dimension.

[smax,detprobsgivesmax,numdetprobs] = calcclassicalbound(n,d,m,corrcoefflist);
darray = zeros(1,n);
marray = zeros(1,n);
probdistsgivesmax = zeros(numdetprobs,(d*m)^n);
for row = 1:numdetprobs
    vars = detprobsgivesmax(row,:);
    probdist = calcprobdist(vars,n,n,darray,marray,n,d,m);
    probdistsgivesmax(row,:) = probdist;
end
reducedmatrix = rref(probdistsgivesmax);
dimension = rank(reducedmatrix);