% A selection of test states for the calcdimandclassicalbound function with the
% expected dimension and maximum classical value achievable, these are found from
% literature. So far the algorithm works for all but the n = 4 tests.

% Test 1: Expected dimension: 8 Expected Bound: 2
% [dimension,classicalbound] = calcdimandclassicalbound(2,2,2,[0 0 0 0 1 1 0 1 -1])

% Test 2: Expected dimension: 15 Expected Bound: 4
% [dimension,classicalbound] = calcdimandclassicalbound(2,2,3,[0 1 1 0 1 -1 -1 1 1 -1 -1 -1 0 1 -1 0])

% Test 3: Expected dimension: 26 Expected Bound: 6
% [dimension,classicalbound] = calcdimandclassicalbound(3,2,2,[+0 +0 +2 +0 +1 -1 +2 -1 -1 +0 +1 -1 +1 -1 -2 -1 -2 +1 +2 -1 -1 -1 -2 +1 -1 +1 +2])

% Test 4: Expected dimension: 63 Expected Bound: 8
% [dimension,classicalbound] = calcdimandclassicalbound(3,2,3,[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 -2 0 0 -2 1 1 0 0 1 -1 0 0 0 0 0 -2 1 1 0 1 4 1 0 1 1 0 0 0 0 0 0 0 1 -1 0 1 1 0 0 -1 0 -1])

% Test 5: Expected dimension: 80 Expected Bound: 9
% [dimension,classicalbound] = calcdimandclassicalbound(4,2,2,[+0 2 1 2 0 0 1 0 -1 2 0 0 0 1 1 0 1 -1 1 0 -1 0 1 -1 -1 -1 0 2 0 0 0 1 1 0 -1 -1 0 1 0 1 -5 2 -1 2 1 0 1 -1 -1 2 1 -1 1 0 1 0 -1 0 1 -1 -1 -1 0 0 1 -1 -1 2 1 -1 1 0 -1 -1 0 -1 1 0 0 0 -2])

% Test 6: Expected dimension: 80 Expected Bound: 10
% [dimension,classicalbound] = calcdimandclassicalbound(4,2,2,[0 3 1 3 0 -1 1 -1 0 3 0 -1 0 0 0 -1 0 -1 1 -1 0 -1 0 -1 0 -1 -1 3 0 -1 0 0 0 -1 0 -1 0 0 0 0 -3 1 0 1 1 -1 0 -1 0 1 1 -1 1 0 1 -1 0 -1 0 -1 0 -1 -1 -1 0 -1 0 1 1 -1 1 0 0 -1 -1 -1 1 0 -1 0 -1])

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
        estimatedtime =  log2((2^(n*m*d))*(((m+1)^n)-1)*((m/(m+1))^2)*d);
        estimatedtimes(i1) = estimatedtime
    end
%%

% Calculate a least squares fit and plot the data and the fit.
p = polyfit(estimatedtimes,calctimes,1)
xfit = linspace(min(estimatedtimes),max(estimatedtimes));
yfit = polyval(p,xfit);
s = scatter(estimatedtimes,calctimes);
s.Marker = '.'; 
xlabel('$\log_{2}(2^{nmd}((m+1)^{n}-1)(\frac{m}{m+1})^{2}d$','Interpreter','Latex','FontSize',15);
ylabel('$\log_{2}(\frac{\Delta{t}}{k})$','Interpreter','Latex','FontSize',15);
hold on
p = plot(xfit,yfit);
p.Color = 'green';
hold off
saveas(gcf,'timings.png')
