function [dimension,smax] = calcdimandclassicalbound(n,d,m,corrcoefflist)
% CALCDIMANDCLASSICALBOUND Calculates the classical bound of a bell
% inequality and it's dimension.
instance = classicalboundcalculator(n,d,m,corrcoefflist);
[smax,detprobsgivesmax] = calcclassicalbound(instance)
darray = zeros(1,n);
marray = zeros(1,n);
numdetprobs = size(detprobsgivesmax,1);
probdistsgivesmax = zeros(numdetprobs,(d*m)^n);
for row = 1:numdetprobs
    vars = detprobsgivesmax(row,:);
    probdist = calcprobdist(vars,n,n,darray,marray,n,d,m);
    probdistsgivesmax(row,:) = probdist;
end
probdistsgivesmax
reducedmatrix = rref(probdistsgivesmax)
dimension = rank(reducedmatrix);