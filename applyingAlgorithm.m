fprintf("Initialising...\n")

% There are two qubits so n = 2, and there are two possible measurement
% operators so m = 2 and there are two possible outcomes for each
% measurement so d = 2.
n = 2;
d = 2;
m = 2;

projs = zeros(2,2,m*d);

%Define Alice and Bobs local measurements, they make the same measurements
phi1 = [sqrt(2)/2;sqrt(2)/2]
phi2 = [sqrt(2)/2;(sqrt(2)/2)*(cos(pi/3)+1i*sin(pi/3))]

%phi1 = [sqrt(2)/2;(sqrt(2)/2)*(cos(pi/8)+1i*sin(pi/8))]
%phi2 = [sqrt(2)/2;(sqrt(2)/2)*(cos(5*pi/8)+1i*sin(5*pi/8))]

%phi1 = [sqrt(2)/2;(sqrt(2)/2)*(cos(7*pi/8)+1i*sin(7*pi/8))]
%phi2 = [sqrt(2)/2;(sqrt(2)/2)*(cos(3*pi/8)+1i*sin(3*pi/8))]

%phi1 = [sqrt(2)/2;(sqrt(2)/2)*(cos(0*pi/8)+1i*sin(0*pi/8))]
%phi2 = [sqrt(2)/2;(sqrt(2)/2)*(cos(pi/2)+1i*sin(pi/2))]

% %Alices measurements
measurements(:,:,1) = phi1*phi1';
measurements(:,:,2) = phi2*phi2';

% %Bobs measurements
measurements(:,:,3) = phi1*phi1';
measurements(:,:,4) = phi2*phi2';

%Calculate the eigenvectors and projectors of Alices measurements
[eigenvectorsMatrix,~] = eig(measurements(:,:,1));
projs(:,:,1) = eigenvectorsMatrix(:,1)*eigenvectorsMatrix(:,1)'
projs(:,:,2) = eigenvectorsMatrix(:,2)*eigenvectorsMatrix(:,2)'
[eigenvectorsMatrix,~] = eig(measurements(:,:,2));
projs(:,:,3) = eigenvectorsMatrix(:,1)*eigenvectorsMatrix(:,1)'
projs(:,:,4) = eigenvectorsMatrix(:,2)*eigenvectorsMatrix(:,2)'

%Calculate the eigenvectors and projectors of Bobs measurements
[eigenvectorsMatrix,~] = eig(measurements(:,:,3));
projs(:,:,5) = eigenvectorsMatrix(:,1)*eigenvectorsMatrix(:,1)'
projs(:,:,6) = eigenvectorsMatrix(:,2)*eigenvectorsMatrix(:,2)'
[eigenvectorsMatrix,~] = eig(measurements(:,:,4));
projs(:,:,7) = eigenvectorsMatrix(:,1)*eigenvectorsMatrix(:,1)'
projs(:,:,8) = eigenvectorsMatrix(:,2)*eigenvectorsMatrix(:,2)'

% Create an array of all the possible tensor products of these projectors
% ensuring that the probabilities in b will be in the right order.
tensorProdProjs = zeros(2^n,2^n,(m*d)^n);
ind = 1;
for i1 = [1,3]
    for i2 = [5,7]
        tensorProdProjs(:,:,ind) = kron(projs(:,:,i1),projs(:,:,i2));
        ind = ind+1;
    end
end
for i1 = [1,3]
    for i2 = [6,8]
        tensorProdProjs(:,:,ind) = kron(projs(:,:,i1),projs(:,:,i2));
        ind = ind+1;
    end
end
for i1 = [2,4]
    for i2 = [5,7]
        tensorProdProjs(:,:,ind) = kron(projs(:,:,i1),projs(:,:,i2));
        ind = ind+1;
    end
end
for i1 = [2,4]
    for i2 = [6,8]
        tensorProdProjs(:,:,ind) = kron(projs(:,:,i1),projs(:,:,i2));
        ind = ind+1;
    end
end

%Create the measurement outcome list for the calculation of A and facet
%dimension of the inequalities
dListElement = (ones(1,m)*d);
for i1 = 1:n
    maxNoMeasOutcomesList(:,i1) = {dListElement'};
end

% Calculate the matrix A
instance = extremalPointLooper(maxNoMeasOutcomesList);
C = calc(instance);
%Transform C into the right form (into A)
A = C;
A(1,:) = C(1,:);
A(2,:) = C(3,:);
A(3,:) = C(9,:);
A(4,:) = C(11,:);
A(5,:) = C(2,:);
A(6,:) = C(4,:);
A(7,:) = C(10,:);
A(8,:) = C(12,:);
A(9,:) = C(5,:);
A(10,:) = C(7,:);
A(11,:) = C(13,:);
A(12,:) = C(15,:);
A(13,:) = C(6,:);
A(14,:) = C(8,:);
A(15,:) = C(14,:);
A(16,:) = C(16,:);

% Define the best possible dimension of the Bell Inequality for the scenario.
bestdimension = ((m*(d-1)+1)^n) - 2;

fprintf("Starting...\n")

% Looping over each state solve Ax=b finding the solution which minimizes 
% the L1 norm of x. Calculate the corresponding Bell Inequality y. 
% Calculate the facet dimension and maximum classical correlation using 
% the facet dimension algorithm   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numPoints = 1; %The number of values for each angle (theta1,theta2,...)
interval = (2*pi/numPoints); %The interval between these values
angleList = linspace(0,2*pi,numPoints);
row = 1; %Initialise the row
totalNumPoints = numPoints^6; %The total number of states to loop over

%Create arrays to store the data
robustnessMatrix = zeros(totalNumPoints,1);
ketStateMatrix = zeros(totalNumPoints,2^d);
inequalityMatrix = zeros(totalNumPoints,(m*d)^n);
dimensionMatrix = zeros(totalNumPoints,1);
bestDimensionMatrix = zeros(totalNumPoints,1);
dimDiffMatrix = zeros(totalNumPoints,1);
classicalBoundMatrix = zeros(totalNumPoints,1);
sMaxMatrix = zeros(totalNumPoints,1);
anglesMatrix = zeros(totalNumPoints,6);

%Looping over all the different values of the parameters
for theta1 = angleList
    for theta2 = angleList
        for theta3 = angleList
            for phi1 = angleList
                for phi2 = angleList
                    for phi3 = angleList
                        
                        %Calculate the complex coefficients
                        c1 = cos(theta1)*cos(theta2)*cos(theta3);
                        c2 = exp(1i*phi1)*(cos(theta2)*cos(theta3)*sin(theta1));
                        c3 = exp(1i*phi2)*(cos(theta3)*sin(theta2));
                        c4 = exp(1i*phi3)*sin(theta3);
                        
                        %Store the angles
                        angles = [theta1,theta2,theta3,phi1,phi2,phi3];
                        ketstate = [c1;c2;c3;c4];
                        
                        %Fix the state for tests
                        ketstate = [1/sqrt(2);0;0;1/sqrt(2)];
                        
                        rho = ketstate*ketstate';
                        numelx = size(tensorProdProjs,3);
                       

                                            
                        %Calculate the vector b
                        b = zeros(numelx,1);
                        for i1 = 1:numelx
                            b(i1) = real(trace(rho*tensorProdProjs(:,:,i1)));
                        end                   

                        % Start timing the calculation.
                        t1 = clock;
    
                        % Apply the CVX optimization algorithm and calculate the robustness and corresponding Bell Inequality y.
                        cvx_begin quiet;
                            variable x(numelx);
                            dual variable y;
                            minimize(norm(x,1));
                            subject to;
                                y: A*x == b;
                        cvx_end;
                        robustness = cvx_optval;

                        %Renormalise these inequalities
                        maxY = max(abs(y.'));
                        yNormalised = y.'/maxY;
                        
                        %Transform the inequality into the right form for
                        %input into the facet dimension algorithm
                        yNormalisedNew(1) = yNormalised(1);
                        yNormalisedNew(5) = yNormalised(2);
                        yNormalisedNew(9) = yNormalised(3);
                        yNormalisedNew(13) = yNormalised(4);
                        yNormalisedNew(2) = yNormalised(5);
                        yNormalisedNew(6) = yNormalised(6);
                        yNormalisedNew(10) = yNormalised(7);
                        yNormalisedNew(14) = yNormalised(8);
                        yNormalisedNew(3) = yNormalised(9);
                        yNormalisedNew(7) = yNormalised(10);
                        yNormalisedNew(11) = yNormalised(11);
                        yNormalisedNew(15) = yNormalised(12);
                        yNormalisedNew(4) = yNormalised(13);
                        yNormalisedNew(8) = yNormalised(14);
                        yNormalisedNew(12) = yNormalised(15);
                        yNormalisedNew(16) = yNormalised(16);
                                        
                        % Apply the facet dimension algorithm to this inequality and calculate its dimension
                        probCoeffList = createProbCoeffList(n,d,m,yNormalisedNew);
                        [dimension,smax] = calcdimandclassicalbound(maxNoMeasOutcomesList,probCoeffList);

                        % Store the data
                        ketStateMatrix(row,:) = ketstate;
                        inequalityMatrix(row,:) = yNormalisedNew;
                        robustnessMatrix(row,:) = robustness;
                        dimensionMatrix(row,:) = dimension;
                        bestDimensionMatrix(row,:) = bestdimension;
                        dimDiffMatrix(row,:) = dimension - bestdimension;
                        sMaxMatrix(row,:) = smax;
                        anglesMatrix(row,:) = angles;

                        T.Properties.RowNames = row;

                        % Stop timing the calculation.
                        t2 = clock;
                        % Calculate the time taken for this calculation.
                        numsecs = etime(t2,t1);
                        % Add this to an array.
                        numsecsarray(row) = numsecs;
                        % Calculate the estimated time left in seconds, from the average time taken per calculation.
                        timeleftsecs = (totalNumPoints-row)*mean(numsecsarray);
                        row = row + 1;
                        % Print out the estimated time left.
                        estimatedtimeleft = datestr(timeleftsecs/(24*60*60), 'DD:HH:MM:SS.FFF')
                        row
                    end
                end
            end
        end
    end
end

%Plot the dimension difference for each of these inequalities
stateParametrisation = 1:numPoints^6;
s2 = scatter(stateParametrisation,dimDiffMatrix,100);
s2.Marker = '.';
xlabel('State','Interpreter','Latex','FontSize',15);
ylabel('$dim - dim_{max}$','Interpreter','Latex','FontSize',15);

minSMax = 999;
maxSMax = 0;
%Calculate how many of these have a dimension of 8
inequalities = inequalityMatrix(dimDiffMatrix==1,:);
numTrivialInequalities = size(inequalities,1)
%Calculate how many of these are unique
uniqueInequalities = unique(round(inequalities,4),'rows');
numUniqueTrivialInequalities = size(uniqueInequalities,1)
%Check how s varies for the 8 dimensional inequalities
minSMax = min(sMaxMatrix(dimDiffMatrix==1))
maxSMax = max(sMaxMatrix(dimDiffMatrix==1))

minSMax = 999;
maxSMax = 0;
%Calculate how many of these have a dimension == maxdimension (7)
inequalities = inequalityMatrix(dimDiffMatrix==0,:);
numTightInequalities = size(inequalities,1)
%Calculate how many of these are unique
uniqueInequalities = unique(round(inequalities,4),'rows');
numUniqueTightInequalities = size(uniqueInequalities,1)
%Check how s varies for these tight Bell Inequalities
minSMax = min(sMaxMatrix(dimDiffMatrix==0))
maxSMax = max(sMaxMatrix(dimDiffMatrix==0))

