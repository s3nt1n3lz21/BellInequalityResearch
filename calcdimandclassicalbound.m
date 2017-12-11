function [dimension,smax] = calcdimandclassicalbound(dlist,coefflist)
% CALCDIMANDCLASSICALBOUND Calculates the classical bound of a bell
% inequality and it's dimension.

% Create an instance of a classicalboundanddimensioncalculator
instance = cbanddimcalc(dlist,coefflist);
% Calculate the dimension and classical bound
[dimension,smax] = calc(instance);