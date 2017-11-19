function [dimension,smax] = calcdimandclassicalbound(n,dlist,corrcoefflist)
% CALCDIMANDCLASSICALBOUND Calculates the classical bound of a bell
% inequality and it's dimension.

% Create an instance of a classicalboundanddimensioncalculator
instance = cbanddimcalc(n,dlist,corrcoefflist);
% Calculate the dimension and classical bound
[dimension,smax] = calc(instance);