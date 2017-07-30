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
% 2^(nmd) * (m+1)^n.
estimatedtimes = zeros(1,numtests);
numtestruns = 2;
timesarray = zeros(numtestruns,numtests);
for i2 = 1:numtestruns
    for i1 = 1:numtests
        coefflist = cell2mat(listofcoefflist(i1));
        scenario = scenarios(i1,:);
        n = scenario(1);
        d = scenario(2);
        m = scenario(3);
        t1 = clock;
        [dimension,classicalbound] = calcdimandclassicalbound(n,d,m,coefflist);
        t2 = clock;
        numsecs = etime(t2,t1);
        timesarray(i2,i1) = log2(1000*numsecs)
        estimatedtime =  n*(m*d+log2(m+1));
        estimatedtimes(i1) = estimatedtime
    %     log2(n*m*d)
    end
end
averagetimes = mean(timesarray,1)
%% 

% Calculate a least squares fit and plot the data and the fit.
p = polyfit(estimatedtimes,averagetimes,1)
xfit = linspace(min(estimatedtimes),max(estimatedtimes));
yfit = polyval(p,xfit);
s = scatter(estimatedtimes,averagetimes);
s.Marker = '.'; 
xlabel('log_{2}(2^{nmd}(m+1)^{n})');
ylabel('log_{2}(\Deltat)');
hold on
p = plot(xfit,yfit);
p.Color = 'green';
hold off
savefig('Timings.fig')

% f = @() myComputeFunction; % handle to function
% timeit(f)
% Look at better ways of timing the algorithm, i.e cputime not matlab time
% etc..
