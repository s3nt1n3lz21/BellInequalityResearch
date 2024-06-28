% A selection of test states for the calcdimandclassicalbound function with the
% expected dimension and maximum classical value achievable, these are found from
% literature. So far the algorithm works for all but the n = 4 tests.

% Test 1: Expected dimension: 7 Expected Bound: 2
% OLD VERSION
% [maxNoMeasOutcomesList,probCoeffList] = convertCorrToProb(2,2,2,[0 0 0 0 1 1 0 1 -1])
% [dimension,classicalbound] = calcdimandclassicalbound(maxNoMeasOutcomesList,probCoeffList)
% NEW VERSION
% maxNoMeasOutcomesList(:,1) = {[2;2]}; maxNoMeasOutcomesList(:,2) = {[2;2]};
% probCoeffList(:,1) = {[0]}; probCoeffList(:,2) = {[0; 0]}; probCoeffList(:,3) = {[0; 0]}; probCoeffList(:,4) = {[0; 0]}; probCoeffList(:,5) = {[1; -1 ;-1 ;1]}; probCoeffList(:,6) = {[1; -1 ;-1; 1]}; probCoeffList(:,7) = {[0; 0]}; probCoeffList(:,8) = {[1; -1; -1; 1]}; probCoeffList(:,9) = {[-1; 1; 1; -1]};
% [dimension,classicalbound] = calcdimandclassicalbound(maxNoMeasOutcomesList,probCoeffList)

% Test 2 : Expected dimension: 14 Expected bound: 0
% maxNoMeasOutcomesList(:,1) = {[2;2;2]}; maxNoMeasOutcomesList(:,2) = {[2;2;2]};
% probCoeffList(:,1) = {[0]};
% probCoeffList(:,2) = {[-2;0]};
% probCoeffList(:,3) = {[-1;0]};
% probCoeffList(:,4) = {[0;0]};
% probCoeffList(:,5) = {[-1;0]};
% probCoeffList(:,6) = {[1;0;0;0]};
% probCoeffList(:,7) = {[1;0;0;0]};
% probCoeffList(:,8) = {[1;0;0;0]};
% probCoeffList(:,9) = {[0;0]};
% probCoeffList(:,10) = {[1;0;0;0]};
% probCoeffList(:,11) = {[1;0;0;0]};
% probCoeffList(:,12) = {[-1;0;0;0]};
% probCoeffList(:,13) = {[0;0]};
% probCoeffList(:,14) = {[1;0;0;0]};
% probCoeffList(:,15) = {[-1;0;0;0]};
% probCoeffList(:,16) = {[0;0;0;0]};

% Test 3: Expected dimension: 25 Expected Bound: 6
% [maxNoMeasOutcomesList,probCoeffList] = convertCorrToProb(3,2,2,[+0 +0 +2 +0 +1 -1 +2 -1 -1 +0 +1 -1 +1 -1 -2 -1 -2 +1 +2 -1 -1 -1 -2 +1 -1 +1 +2])
% [dimension,classicalbound] = calcdimandclassicalbound(maxNoMeasOutcomesList,probCoeffList)

% Test 4: Expected dimension: 62 Expected Bound: 8
% [maxNoMeasOutcomesList,probCoeffList] = convertCorrToProb(3,2,3,[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 -2 0 0 -2 1 1 0 0 1 -1 0 0 0 0 0 -2 1 1 0 1 4 1 0 1 1 0 0 0 0 0 0 0 1 -1 0 1 1 0 0 -1 0 -1])
% [dimension,classicalbound] = calcdimandclassicalbound(maxNoMeasOutcomesList,probCoeffList)

% Test 5: Expected dimension: 79 Expected Bound: 2
% [maxNoMeasOutcomesList,probCoeffList] = convertCorrToProb(4,2,2,[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 -1])
% [dimension,classicalbound] = calcdimandclassicalbound(maxNoMeasOutcomesList,probCoeffList)

% Test 6: Expected dimension: 241 Expected Bound: 2
% [maxNoMeasOutcomesList,probCoeffList] = convertCorrToProb(5,2,2,[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1])
% [dimension,classicalbound] = calcdimandclassicalbound(maxNoMeasOutcomesList,probCoeffList)

% Test 7: Expected dimension: 10 Expected Bound: 0 
% maxNoMeasOutcomesList(:,1) = {[2;2]}; maxNoMeasOutcomesList(:,2) = {[2;2;2]};
% probCoeffList(:,1) = {[0]};
% probCoeffList(:,2) = {[0;0]};
% probCoeffList(:,3) = {[-1;0]};
% probCoeffList(:,4) = {[0;0]};
% probCoeffList(:,5) = {[-1;0]};
% probCoeffList(:,6) = {[0;0;0;0]};
% probCoeffList(:,7) = {[1;0;0;0]};
% probCoeffList(:,8) = {[1;0;0;0]};
% probCoeffList(:,9) = {[0;0]};
% probCoeffList(:,10) = {[0;0;0;0]};
% probCoeffList(:,11) = {[1;0;0;0]};
% probCoeffList(:,12) = {[-1;0;0;0]};

% Test ?: Expected dimension: 23 Expected Bound: 0 equation 39 doesnt work?
% no terms for third outcome? reading incorrectly? talking about
% correlators.
% maxNoMeasOutcomesList(:,1) = {[2;2]}; maxNoMeasOutcomesList(:,2) = {[2;2]};
% probCoeffList(:,1) = {[0]};
% probCoeffList(:,2) = {[-1;-1;0]};
% probCoeffList(:,3) = {[0;0;0]};
% probCoeffList(:,4) = {[-1;-1;0]};
% probCoeffList(:,5) = {[1;1;0;1;0;0;0;0;0]};
% probCoeffList(:,6) = {[0;1;0;1;1;0;0;0;0]};
% probCoeffList(:,7) = {[0;0;0]};
% probCoeffList(:,8) = {[0;1;0;1;1;0;0;0;0]};
% probCoeffList(:,9) = {[0;-1;0;-1;-1;0;0;0;0]};

% probCoeffList(:,1) = {[0]};
% probCoeffList(:,2) = {[-1;-1]};
% probCoeffList(:,3) = {[0;0]};
% probCoeffList(:,4) = {[-1;-1]};
% probCoeffList(:,5) = {[1;1;1;0]};
% probCoeffList(:,6) = {[0;1;1;1]};
% probCoeffList(:,7) = {[0;0]};
% probCoeffList(:,8) = {[0;1;1;1]};
% probCoeffList(:,9) = {[0;-1;-1;-1]};

% Written down correctly?
% equation 14 in "relevant multi-setting..." dimension:23 bound:2
% maxNoMeasOutcomesList(:,1) = {[3;3]}; maxNoMeasOutcomesList(:,2) = {[3;3]};
% probCoeffList(:,1) = {[0]};
% probCoeffList(:,2) = {[0;0;0]};
% probCoeffList(:,3) = {[0;0;0]};
% probCoeffList(:,4) = {[0;0;0]};
% probCoeffList(:,5) = {[1;0;-1;-1;1;0;0;-1;1]};
% probCoeffList(:,6) = {[1;-1;0;0;1;-1;-1;0;1]};
% probCoeffList(:,7) = {[0;0;0]};
% probCoeffList(:,8) = {[-1;0;1;1;-1;0;0;1;-1]};
% probCoeffList(:,9) = {[1;0;-1;-1;1;0;0;-1;1]};

% Test 2 old: Expected dimension: 14 Expected Bound: 4 reference?
% OLD VERSION
% [maxNoMeasOutcomesList,probCoeffList] = convertCorrToProb(2,2,3,[0 1 1 0 1 -1 -1 1 1 -1 -1 -1 0 1 -1 0])
% [dimension,classicalbound] = calcdimandclassicalbound(maxNoMeasOutcomesList,probCoeffList)
% NEW VERSION
% maxNoMeasOutcomesList(:,1) = {[2;2;2]}; maxNoMeasOutcomesList(:,2) = {[2;2;2]};
% probCoeffList(:,1) = {[0]}; probCoeffList(:,2) = {[1 ;-1]}; probCoeffList(:,3) = {[1; -1]}; probCoeffList(:,4) = {[0; 0]}; probCoeffList(:,5) = {[1; -1]}; probCoeffList(:,6) = {[-1; 1; 1; -1]}; probCoeffList(:,7) = {[-1; 1 ;1 ;-1]}; probCoeffList(:,8) = {[1; -1; -1; 1]}; probCoeffList(:,9) = {[1 ;-1]}; probCoeffList(:,10) = {[-1; 1 ;1 ;-1]}; probCoeffList(:,11) = {[-1; 1; 1; -1]}; probCoeffList(:,12) = {[-1; 1; 1; -1]}; probCoeffList(:,13) = {[0; 0]}; probCoeffList(:,14) = {[1; -1; -1; 1]}; probCoeffList(:,15) = {[-1; 1; 1; -1]}; probCoeffList(:,16) = {[0; 0; 0; 0]};
% [dimension,classicalbound] = calcdimandclassicalbound(maxNoMeasOutcomesList,probCoeffList)

% INCORRECT TEST
% Test ?: Expected dimension: 80 Expected Bound: 9
% [maxNoMeasOutcomesList,probCoeffList] = convertCorrToProb(4,2,2,[+0 2 1 2 0 0 1 0 -1 2 0 0 0 1 1 0 1 -1 1 0 -1 0 1 -1 -1 -1 0 2 0 0 0 1 1 0 0 -1 0 1 0 1 -5 2 -1 2 1 0 1 -1 -1 2 1 -1 1 0 1 0 -1 0 1 -1 -1 -1 0 0 1 -1 -1 2 1 -1 1 0 -1 -1 0 -1 1 0 0 0 -2])
% [dimension,classicalbound] = calcdimandclassicalbound(maxNoMeasOutcomesList,probCoeffList)

% INCORRECT TEST
% Test ?: Expected dimension: 80 Expected Bound: 10
% [maxNoMeasOutcomesList,probCoeffList] = convertCorrToProb(4,2,2,[0 3 1 3 0 -1 1 -1 0 3 0 -1 0 0 0 -1 0 -1 1 -1 0 -1 0 -1 0 -1 -1 3 0 -1 0 0 0 -1 0 -1 0 0 0 0 -3 1 0 1 1 -1 0 -1 0 1 1 -1 1 0 1 -1 0 -1 0 -1 0 -1 -1 -1 0 -1 0 1 1 -1 1 0 0 -1 -1 -1 1 0 -1 0 -1])
% [dimension,classicalbound] = calcdimandclassicalbound(maxNoMeasOutcomesList,probCoeffList)

% Code to test the timings of the algorithm using these test states, here
% only the test states that the algorithm succeeded producing the correct
% results for will be used.

% Create a cell array to hold the coefficient lists for each test.
listofcoefflist = {[0 0 0 0 1 1 0 1 -1];
                   [0 1 1 0 1 -1 -1 1 1 -1 -1 -1 0 1 -1 0];
                   [+0 +0 +2 +0 +1 -1 +2 -1 -1 +0 +1 -1 +1 -1 -2 -1 -2 +1 +2 -1 -1 -1 -2 +1 -1 +1 +2];
                   [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 -2 0 0 -2 1 1 0 0 1 -1 0 0 0 0 0 -2 1 1 0 1 4 1 0 1 1 0 0 0 0 0 0 0 1 -1 0 1 1 0 0 -1 0 -1]};
% Create an array to hold the corresponding scenarios.
scenarios = [[2,2,2];[2,2,3];[3,2,2];[3,2,3]];
numtests = size(scenarios,1);

% Create an array to hold the estimated time taken for each test. deltat ~
% 2^(nmd) * ((m+1)^n-1).
estimatedtimes = zeros(1,numtests);
calctimes = zeros(1,numtests);
    for i1 = 1:numtests
        coefflist = cell2mat(listofcoefflist(i1));
        scenario = scenarios(i1,:);
        n = scenario(1);
        d = scenario(2);
        m = scenario(3);
        f = @() calcdimandclassicalbound(n,d,m,coefflist);
        t = timeit(f,2)
        nnzterms = nnz(coefflist);
        fractionfull = nnzterms/(((m+1)^n)-1);
        numsecs = t/fractionfull;
        calctimes(i1) = log2(1000*numsecs)
        estimatedtime =  log2((d^(n*m))*(((m+1)^n)-1)*((m/(m+1))^2)*d);
        estimatedtimes(i1) = estimatedtime
    end
%%

% Calculate a least squares fit and plot the data and the fit.
[p,s] = polyfit(estimatedtimes,calctimes,1)
xfit = linspace(min(estimatedtimes),max(estimatedtimes));
[yfit,yfiterr] = polyval(p,xfit,s);
s = scatter(estimatedtimes,calctimes);
s.Marker = '.'; 
xlabel('$\log_{2}\left(d^{nm}\left(\left(m+1\right)^{n}-1\right)\left(\frac{m}{m+1}\right)^{2}d\right)$','Interpreter','Latex','FontSize',15);
ylabel('$\log_{2}\left(\frac{\Delta{t}}{k}\right)$','Interpreter','Latex','FontSize',15);
hold on
plot(xfit,yfit,['--',[],'g'])
hold off
saveas(gcf,'timings.png')
