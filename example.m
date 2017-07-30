% An example of the application of the algorithm to some states. In this
% file we calculate the robustness and Bell inequality of each state and
% then apply this algorithm to these Bell inequalities. The dimension is compared to the best
% possible dimension and an estimation of the quantum bound is also calculated
% from the maximum eigenvalue of the bell operator. The corresponding
% eigenvector is also calculated but this is left commented out.

% The states given in the file are 3 qubit states and each row represents
% one state with each element being the coefficient in the standard basis.
% To calculate the robustness and the corresponding bell inequality we have to
% find the decomposition of the density matrix rho, into the basis of tensor products
% of projector states.

fprintf("Initialising...\n")

% First define the pauli operators.
paulis = cat(3,[1,0;0,1],[0,1;1,0],[0,-1i;1i,0],[1,0;0,-1]);

% Create the six projection operators from these.
for i1 = 1:3
    sproj(:,:,i1) = 0.5*(eye(2)+paulis(:,:,i1+1));
    sproj(:,:,i1+3) = 0.5*(eye(2)-paulis(:,:,i1+1));
end

% There are three qubits so n = 3, and there are three pauli measurement
% operators so m = 3 and there are two possible outcomes for each
% measurement so d = 2.
n = 3;
d = 2;
m = 3;

% Create an array of all the possible tensor products of these projectors
% and their decomposition into the basis.
tensorprojs = zeros(2^n,2^n,6^n);
vectensorprojs = zeros((m+1)^n,6^n);
ind = 1;
for i1 = 1:6
    for i2 = 1:6
        for i3 = 1:6
            tensorprojs(:,:,ind) = kron(sproj(:,:,i1),kron(sproj(:,:,i2),sproj(:,:,i3)));
            vectensorprojs(:,ind) = rhodecomp(tensorprojs(:,:,ind));
            ind = ind+1;
        end
    end
end

% Define the best possible dimension for the scenario.
bestdimension = ((m*(d-1)+1)^n) - 1;
% Get the states from file.
states = importdata('Three_qubit_hierarchy.mat');
numrows = size(states,1);

% Create arrays to hold all of the data.
statematrix = zeros(numrows,8);
robustnessmatrix = zeros(numrows,1);
inequalitymatrix = zeros(numrows,(m+1)^n);
dimensionmatrix = zeros(numrows,1);
bestdimensionmatrix = zeros(numrows,1);
dimdiffmatrix = zeros(numrows,1);
classicalboundmatrix = zeros(numrows,1);
quantumboundmatrix = zeros(numrows,1);
bounddiffmatrix = zeros(numrows,1);
% beststatematrix = zeros(numrows,8);

fprintf("Starting...\n")

% Looping over each state calculate the decomposition of the state and
% solve Ax=b finding the solution which minimizes the L1 norm of x. The
% minimum L1 norm of x is the robustness of the state. Due to duality, from this we can
% also calculate the corresponding Bell Inequality y. This is then plugged
% into the algorithm to calculate its dimension and classical bound.

for row = 1:numrows
    ketstate = states(row,:).';
    b = rhodecomp(ketstate * ketstate').';
    numelx = size(tensorprojs,3);
    A = vectensorprojs;

    % Start timing the calculation.
    t1 = clock;
    
    % Calculate the robustness and corresponding Bell Inequality y.
    cvx_begin quiet; 
        variable x(numelx);
        dual variable y;
        minimize( norm( x, 1 ) );
        subject to;
          y: A * x == b;
    cvx_end;
    robustness=cvx_optval;
    
    % Apply the algorithm.
    [dimension,smax] = calcdimandclassicalbound(n,d,m,y.');
    
    % Store the data.
    robustnessmatrix(row,:) = norm(x,1);
    inequalitymatrix(row,:) = y.';
    dimensionmatrix(row,:) = dimension;
    bestdimensionmatrix(row,:) = bestdimension;
    dimdiffmatrix(row,:) = dimension - bestdimension;
    classicalboundmatrix(row,:) = smax;
    
    % Calculate the Bell Operator b to calculate the estimation of the
    % quantum bound.
    b = 0;
    for index = 1:length(inequalitymatrix(row,:))
        settings = dec2base(index-1,m+1,n)-'0';
        numset = length(settings);
        oparray = zeros(2,2,numset);
        for party = 1: numset
            settingindex = settings(party);
            operator = paulis(:,:,settingindex+1);
            oparray(:,:,party) = operator;
        end
        tensprod = kron(oparray(:,:,1),oparray(:,:,2));
        b = b + inequalitymatrix(row,index)*kron(tensprod,oparray(:,:,3));
    end
    
    % Calculate the eigenvalue.
    [matrixeigvec,matrixeig] = eig(b);
    
    eigenvalues = zeros(1,length(matrixeig));
    for i1 = 1:length(matrixeig)
    eigenvalues(i1) = matrixeig(i1,i1);
    end
    quantum = max(eigenvalues);
%     index = find(eigenvalues == quantum);
%     beststate = matrixeigvec(:,index);
   
    quantumboundmatrix(row,:) = quantum;
%     beststatematrix(row,:) = beststate.';       
    
    T.Properties.RowNames = row;
    
    % End timing the calculation.
    t2 = clock;
    % Calculate the time taken for this calculation.
    numsecs = etime(t2,t1);
    % Add this to an array.
    numsecsarray(row) = numsecs;
    % Calculate the estimated time left in seconds, from the average time taken per calculation.
    timeleftsecs = (numrows-row)*mean(numsecsarray);
    % Print out the estimated time left.
    estimatedtimeleft = datestr(timeleftsecs/(24*60*60), 'DD:HH:MM:SS.FFF')
    
end

% Put all the data into a table.
data = table(states,robustnessmatrix,inequalitymatrix,dimensionmatrix,bestdimensionmatrix,dimdiffmatrix,classicalboundmatrix,quantumboundmatrix,bounddiffmatrix);
