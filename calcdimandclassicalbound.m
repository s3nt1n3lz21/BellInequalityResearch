function [dimension,smax] = calcdimandclassicalbound(n,d,m,corrcoefflist)
% CALCDIMANDCLASSICALBOUND Calculates the classical bound of a bell
% inequality and it's dimension.

% Create an instance of a classicalboundanddimensioncalculator
instance = cbanddimcalc(n,d,m,corrcoefflist);
% Calculate the dimension and classical bound
[dimension,smax] = calc(instance);