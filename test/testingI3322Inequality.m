% In this algorithm we apply the algorithm to the Bell Inequalities found by using some quantum states and measurement projectors. 
% In this file we calculate the robustness and the Bell inequality of each state and
% then apply this algorithm to these Bell inequalities.

% Calculate the matrix A for CHSH by defining a class to loop over all of
% the extremal points and calculate their vectors. The matrix A has columns
% that are the behaviours of these extremal points.
fprintf("Initialising...\n")

% There are two qubits so n = 2, and there are two possible measurement
% operators so m = 2 and there are two possible outcomes for each
% measurement so d = 2.
n = 2;
d = 2;
m = 2;

% Define the projection operators for CHSH. There are d projectors for the
% different outcomes of each measurement. If there are m different types of
% measurement each with a different basis then there are m*d projectors in
% total.
ket0 = [1; 0];
ket1 = [0; 1];
ketPlus = (1/sqrt(2))*(ket0+ket1);
ketMinus = (1/sqrt(2))*(ket0-ket1);
projs = zeros(2,2,m*d);

%Define the sigma measurement operators
sigmaX = ket0*ket1'+ket1*ket0';
sigmaY = [0,-1i;1i,0]
sigmaZ = ket0*ket0'-ket1*ket1';

% %Alices measurements
% projs(:,:,1) = ket0*ket0';
% projs(:,:,2) = ket1*ket1';
% 
% %Bobs measurements
% projs(:,:,3) = ketPlus*ketPlus';
% projs(:,:,4) = ketMinus*ketMinus';

%Define Alices measurements
% measurements(:,:,1) = sigmaZ;
% measurements(:,:,2) = sigmaX;

%Define Bobs measurements
% measurements(:,:,3) = (1/sqrt(2))*(sigmaZ+sigmaX);
% measurements(:,:,4) = (1/sqrt(2))*(sigmaZ-sigmaX);

%Define Alice and Bobs local measurements, they make the same measurements
%phi1 = [sqrt(2)/2;sqrt(2)/2]
%phi2 = [sqrt(2)/2;(sqrt(2)/2)*(cos(pi/3)+1i*sin(pi/3))]

phi1 = [sqrt(2)/2;(sqrt(2)/2)*(cos(pi/8)+1i*sin(pi/8))]
phi2 = [sqrt(2)/2;(sqrt(2)/2)*(cos(5*pi/8)+1i*sin(5*pi/8))]

%phi1 = [sqrt(2)/2;(sqrt(2)/2)*(cos(7*pi/8)+1i*sin(7*pi/8))]
%phi2 = [sqrt(2)/2;(sqrt(2)/2)*(cos(3*pi/8)+1i*sin(3*pi/8))]

%phi1 = [sqrt(2)/2;(sqrt(2)/2)*(cos(0*pi/8)+1i*sin(0*pi/8))]
%phi2 = [sqrt(2)/2;(sqrt(2)/2)*(cos(pi/2)+1i*sin(pi/2))]

%phi1 = 

% %Alices measurements
measurements(:,:,1,1) = phi1*phi1';
measurements(:,:,1,2) = phi2*phi2';

% %Bobs measurements
measurements(:,:,2,1) = phi1*phi1';
measurements(:,:,2,2) = phi2*phi2';

% % %Alices measurements
% measurements(:,:,1,1) = sigmaX;
% measurements(:,:,1,2) = sigmaY;
% measurements(:,:,1,3) = sigmaZ;
% 
% % %Bobs measurements
% measurements(:,:,2,1) = sigmaX;
% measurements(:,:,2,2) = sigmaY;
% measurements(:,:,2,3) = sigmaZ;

%Calculate the eigenvectors and projectors for each parties' measurements
for iN = 1:n
    for iM = 1:m
        [eigenvectorsMatrix,~] = eig(measurements(:,:,1,1));
        for iD = 1:d
            projs(:,:,iN,iD,iM) = eigenvectorsMatrix(:,iD)*eigenvectorsMatrix(:,iD)'
        end
    end
end

% %Calculate the eigenvectors and projectors of Alices measurements
% [eigenvectorsMatrix,~] = eig(measurements(:,:,1,1));
% projs(:,:,1,1,1) = eigenvectorsMatrix(:,1)*eigenvectorsMatrix(:,1)'%*measurements(:,:,1); %|\psi><\psi|A = \lambda|\psi><\psi|
% projs(:,:,1,2,1) = eigenvectorsMatrix(:,2)*eigenvectorsMatrix(:,2)'%*measurements(:,:,1);
% [eigenvectorsMatrix,~] = eig(measurements(:,:,1,2));
% projs(:,:,1,1,2) = eigenvectorsMatrix(:,1)*eigenvectorsMatrix(:,1)'%*measurements(:,:,2);
% projs(:,:,1,2,2) = eigenvectorsMatrix(:,2)*eigenvectorsMatrix(:,2)'%*measurements(:,:,2);
% [eigenvectorsMatrix,~] = eig(measurements(:,:,1,3));
% projs(:,:,1,1,3) = eigenvectorsMatrix(:,1)*eigenvectorsMatrix(:,1)'%*measurements(:,:,2);
% projs(:,:,1,2,3) = eigenvectorsMatrix(:,2)*eigenvectorsMatrix(:,2)'%*measurements(:,:,2);
% 
% %Calculate the eigenvectors and projectors of Bobs measurements
% [eigenvectorsMatrix,~] = eig(measurements(:,:,2,1));
% projs(:,:,2,1,1) = eigenvectorsMatrix(:,1)*eigenvectorsMatrix(:,1)'%*measurements(:,:,3);
% projs(:,:,2,2,1) = eigenvectorsMatrix(:,2)*eigenvectorsMatrix(:,2)'%*measurements(:,:,3);
% [eigenvectorsMatrix,~] = eig(measurements(:,:,2,2));
% projs(:,:,2,1,2) = eigenvectorsMatrix(:,1)*eigenvectorsMatrix(:,1)'%*measurements(:,:,4);
% projs(:,:,2,2,2) = eigenvectorsMatrix(:,2)*eigenvectorsMatrix(:,2)'%*measurements(:,:,4);
% [eigenvectorsMatrix,~] = eig(measurements(:,:,2,3));
% projs(:,:,2,1,3) = eigenvectorsMatrix(:,1)*eigenvectorsMatrix(:,1)'%*measurements(:,:,4);
% projs(:,:,2,2,3) = eigenvectorsMatrix(:,2)*eigenvectorsMatrix(:,2)'%*measurements(:,:,4);

% Create an array of all the possible tensor products of these projectors
% ensuring that the probabilities in b will be in the right order.

D1 = d
D2 = d
M1 = m
M2 = m

tensorProdProjs = zeros(2^n,2^n,(m*d)^n);

ind = 1;
for d1 = 1:D1
    for d2 = 1:D2
        for m1 = 1:M1
            for m2 = 1:M2
                x = [d1,d2,m1,m2];
                display(x)
                tensorProdProjs(:,:,ind) = kron(projs(:,:,1,m1,d1),projs(:,:,2,m2,d2));
                ind = ind + 1;
            end
        end
    end
end

% tensorProdProjs = zeros(2^n,2^n,(m*d)^n);
% ind = 1;
% for i1 = [1,3]
%     for i2 = [5,7]
%         tensorProdProjs(:,:,ind) = kron(projs(:,:,i1),projs(:,:,i2));
%         ind = ind+1;
%     end
% end
% for i1 = [1,3]
%     for i2 = [6,8]
%         tensorProdProjs(:,:,ind) = kron(projs(:,:,i1),projs(:,:,i2));
%         ind = ind+1;
%     end
% end
% for i1 = [2,4]
%     for i2 = [5,7]
%         tensorProdProjs(:,:,ind) = kron(projs(:,:,i1),projs(:,:,i2));
%         ind = ind+1;
%     end
% end
% for i1 = [2,4]
%     for i2 = [6,8]
%         tensorProdProjs(:,:,ind) = kron(projs(:,:,i1),projs(:,:,i2));
%         ind = ind+1;
%     end
% end

%Put the numbers n,m,d into the correct form for input to the
%extremalPointLooper class.
dListElement = (ones(1,m)*d);
for i1 = 1:n
    maxNoMeasOutcomesList(:,i1) = {dListElement'};
end

% Create an instance of a extremalPointLooper class
instance = extremalPointLooper(maxNoMeasOutcomesList);
% Calculate the matrix A
C = calc(instance);
A = C;
% A(3,:) = C(5,:);
% A(4,:) = C(6,:);
% A(5,:) = C(3,:);
% A(6,:) = C(4,:);
% A(11,:) = C(13,:);
% A(12,:) = C(14,:);
% A(13,:) = C(11,:);
% A(14,:) = C(12,:);
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
% Get the state
states = zeros(1,4);
states(1,:) = [0, 1/sqrt(2), -1/sqrt(2), 0];
%states(2,:) = [3/8, 1/8, 1/8, 3/8];
numrows = size(states,1);

% Create arrays to hold all of the data.
%stateMatrix = zeros(numrows,8);

fprintf("Starting...\n")

% Looping over each state calculate the decomposition of the state and
% solve Ax=b finding the solution which minimizes the L1 norm of x. The
% minimum L1 norm of x is the robustness of the state. Due to duality, from this we can
% also calculate the corresponding Bell Inequality y. This is then plugged
% into the algorithm to calculate its dimension and classical bound.    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numPoints = 1;
interval = (2*pi/numPoints);
%angleList = [i*interval for i in range(1,numPoints)]
angleList = linspace(0,2*pi,1);
theta1List = linspace(0,2*pi,numPoints);
pointCounter = 1;
totalNumPoints = numPoints^6;

robustnessMatrix = zeros(totalNumPoints,1);
ketStateMatrix = zeros(totalNumPoints,2^d);
inequalityMatrix = zeros(totalNumPoints,(m*d)^n);
dimensionMatrix = zeros(totalNumPoints,1);
bestDimensionMatrix = zeros(totalNumPoints,1);
dimDiffMatrix = zeros(totalNumPoints,1);
classicalBoundMatrix = zeros(totalNumPoints,1);

%Looping over all the different values of the parameters
for theta1 = theta1List
    for theta2 = theta1List
        for theta3 = theta1List
            for phi1 = theta1List
                for phi2 = theta1List
                    for phi3 = theta1List
                        %theta1;
                        %theta2;
                        %theta3;
                        %phi1;
                        %phi2;
                        %phi3;
                        c1 = cos(theta1)*cos(theta2)*cos(theta3);
                        c2 = exp(1i*phi1)*(cos(theta2)*cos(theta3)*sin(theta1));
                        c3 = exp(1i*phi2)*(cos(theta3)*sin(theta2));
                        c4 = exp(1i*phi3)*sin(theta3);
                        ketstate = [c1;c2;c3;c4];
                        numelx = size(tensorProdProjs,3);
                        b = zeros(numelx,1);
                        ketstate = [1/sqrt(2);0;0;1/sqrt(2)];
                        %ketstate = [cos(pi/8),0,0,sin(pi/8)];
                        rho = ketstate*ketstate';
                        %trace(rho);
                        %tensorProdProjs(:,:,i1)
                        
                        %U = [(cos(pi/16)-1i*sin(pi/16)),0;0,(cos(pi/16)+1i*sin(pi/16))]
                        %rotation = kron(U,U)
                        %rho = ketstate*ketstate'
                        %rho = rotation*rho
                        
                        %Calculate the vector b
                        for i1 = 1:numelx
                            %tensorProdProjs(:,:,i1)
                            %trace(tensorProdProjs(:,:,i1))
                            b(i1) = real(trace(rho*tensorProdProjs(:,:,i1)));
                        end
                        
                        %b = 0.15*b
                        %numelx = 6
                        %b = A(:,8)
                        p = (sqrt(2)+2)/8;
                        %bTest = [p;1/2-p;1/2-p;1/2-p;1/2-p;p;p;p;1/2-p;p;p;p;p;1/2-p;1/2-p;1/2-p];
                        %b = bTest;
                        %A = ATest2;
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

                        %xx
                        %sum(abs(b))
                        % Apply my algorithm to this inequality and
                        % calculate its dimension
                        
                        maxY = max(abs(y.'));
                        yNormalised = y.'/maxY;
                        
                        %yNormalised = fliplr(yNormalisedOrig)
                        %YNormalisedNew = yNormalised;
                        
                        %yNormalised = y.';
                        
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
                        
                        probCoeffList = createProbCoeffList(n,d,m,yNormalisedNew);
                        [dimension,smax] = calcdimandclassicalbound(maxNoMeasOutcomesList,probCoeffList);

                        % Store the data.
                        ketStateMatrix(pointCounter,:) = ketstate;
                        ketstate;
                        pointCounter;
                        ketStateMatrix;
                        inequalityMatrix(pointCounter,:) = yNormalisedNew;
                        robustnessMatrix(pointCounter,:) = robustness;
                        dimensionMatrix(pointCounter,:) = dimension;
                        bestDimensionMatrix(pointCounter,:) = bestdimension;
                        dimDiffMatrix(pointCounter,:) = dimension - bestdimension;   

                        T.Properties.RowNames = pointCounter;

                        % End timing the calculation.
                        t2 = clock;
                        % Calculate the time taken for this calculation.
                        numsecs = etime(t2,t1);
                        % Add this to an array.
                        numsecsarray(pointCounter) = numsecs;
                        % Calculate the estimated time left in seconds, from the average time taken per calculation.
                        timeleftsecs = (totalNumPoints-pointCounter)*mean(numsecsarray);
                        pointCounter = pointCounter + 1;
                        % Print out the estimated time left.
                        estimatedtimeleft = datestr(timeleftsecs/(24*60*60), 'DD:HH:MM:SS.FFF')
                    end
                end
            end
        end
    end
end

% Put all the data into a table.
% bounddiffmatrix = quantumboundmatrix-classicalBoundMatrix;

data = table(robustnessMatrix,inequalityMatrix,dimensionMatrix,bestDimensionMatrix,dimDiffMatrix); % quantumboundmatrix,bounddiffmatrix

%xfit = bounddiffmatrix;
%yfit = dimDiffMatrix;
dimDiffMatrixLessPoly = dimDiffMatrix(dimDiffMatrix<1);
stateParametrisation = 1:size(dimDiffMatrixLessPoly);
s2 = scatter(stateParametrisation,dimDiffMatrixLessPoly,100);
s2.Marker = '.';
xlabel('State','Interpreter','Latex','FontSize',15);
ylabel('$dim - dim_{max}$','Interpreter','Latex','FontSize',15);
inequalities = inequalityMatrix(dimDiffMatrix<1,:);
uniqueInequalities = unique(round(inequalities,4),'rows');
%dimDiffMatrixLess8 = dimDiffMatrix(dimDiffMatrix<1);
%ketStateMatrixUnique = unique(round(ketStateMatrix,4),'rows')
%uniqueInequalityIndices = find()
%inequalityMatrixUnique = unique(round(inequalityMatrix,4),'rows')

%Y = unique(round(inequalityMatrix,2),'rows')
%tightInequalitiesIndices = find(~dimDiffMatrix)
%tightInequalities = inequalityMatrix((tightInequalitiesIndices),:)
%saveas(gcf,'dimvbound.png')

%Testing
%p = (sqrt(2)+2)/8
%bTest = [p;1/2-p;1/2-p;1/2-p;1/2-p;p;p;p;1/2-p;p;p;p;p;1/2-p;1/2-p;1/2-p]
%yNormalised*b =2rrot2 which is right because b is the optimal measurements
%on a maximially entangled state
%chsh = [1 -1 -1 1 1 -1 -1 1 1 -1 -1 1 -1 1 1 -1];
%chsh*b does not produce 2root2 suggesting chsh and b are in different
%orders.
%probCoeffList = createProbCoeffList(n,d,m,chsh);
%%[dimension,smax] = calcdimandclassicalbound(maxNoMeasOutcomesList,probCoeffList);

%y = uniqueInequalities(1,:)
%yNew = y/max(abs(y))
%probCoeffList = createProbCoeffList(n,d,m,yNew)
%[dimension,smax] = calcdimandclassicalbound(maxNoMeasOutcomesList,probCoeffList)


%bellOperator = zeros(2^n,2^n)
%for i1 = 1:size(yNormalised,2)
%    bellOperator = bellOperator + yNormalised(i1)*tensorProdProjs(:,:,i1)
%end
%s = trace(bellOperator*rho)