function [dimension,smax] = calcDimAndClassicalBound(maxNoMeasOutcomesList,probCoeffList)
% CALCDIMANDCLASSICALBOUND Calculates the classical bound of a bell
% inequality and it's dimension.

% Create an instance of a classicalboundanddimensioncalculator
instance = classicalBoundAndDimCalculator(maxNoMeasOutcomesList,probCoeffList);
% Calculate the dimension and classical bound
[dimension,smax] = calc(instance);