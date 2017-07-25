% Here I calculate the robustness and other quantities for a variety of
% states in the file "three_qubit_hierarchary". I calculate the classical
% bound and dimension of the bell inequality found by calculating the
% robustness. I also calculate the difference of the dimension of the bell
% inequality from the dimension of the best possible bell inequality and an
% estimation of the quantum bound using the eigenvalue of the bell
% operator. The corresponding eigenvector represents a better quantum state
% that can be used to test the bell inequality.

% A selection of test states for the calcdimandclassicalbound function with the
% expected dimension and maximum classical value achievable found from
% literature.

% 8,2 (2,2,2,[0 0 0 0 1 1 0 1 -1])
% 15,4 (2,2,3,[0 1 1 0 1 -1 -1 1 1 -1 -1 -1 0 1 -1 0])
% 26,6 (3,2,2,[+0 +0 +2 +0 +1 -1 +2 -1 -1 +0 +1 -1 +1 -1 -2 -1 -2 +1 +2 -1 -1 -1 -2 +1 -1 +1 +2])
% 80,9 (4,2,2,[+0 2 1 2 0 0 1 0 -1 2 0 0 0 1 1 0 1 -1 1 0 -1 0 1 -1 -1 -1 0 2 0 0 0 1 1 0 1 -1 0 1 -1 1 -5 2 -1 2 1 0 1 -1 -1 2 1 -1 1 0 1 0 -1 0 1 -1 -1 -1 0 0 1 -1 -1 2 1 -1 1 0 -1 -1 0 -1 1 0 0 0 -2])
% 80,10 (4,2,2,[0 3 1 3 0 -1 1 -1 0 3 0 -1 0 0 0 -1 0 -1 1 -1 0 -1 0 -1 0 -1 -1 3 0 -1 0 0 0 -1 0 -1 0 0 0 0 -3 1 0 1 1 -1 0 -1 0 1 1 -1 1 0 1 -1 0 -1 0 -1 0 -1 -1 -1 0 -1 0 1 1 - 1 0 0 -1 -1 -1 1 0 -1 0 -1])
% 63,8 (3,2,3,[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 -2 0 0 -2 1 1 0 0 1 -1 0 0 0 0 0 -2 1 1 0 1 4 1 0 1 1 0 0 0 0 0 0 0 1 -1 0 1 1 0 0 -1 0 -1

fprintf("Initialising...\n")

% define the pauli operators
paulis = cat(3,[1,0;0,1],[0,1;1,0],[0,-1i;1i,0],[1,0;0,-1]);

% create the six projection operators from these
for i1 = 1:3
    sproj(:,:,i1) = 0.5*(eye(2)+paulis(:,:,i1+1));
    sproj(:,:,i1+3) = 0.5*(eye(2)-paulis(:,:,i1+1));
end

% create an array of all the possible tensor products of these projectors
% and their decomposition into the basis.

n = 3;
d = 2;
m = 3;
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

bestdimension = 2*(d-1)*m+((d-1)^2)*m^2;
states = importdata('Three_qubit_hierarchy.mat');
numrows = size(states,1);

statematrix = zeros(numrows,8);
robustnessmatrix = zeros(numrows,1);
inequalitymatrix = zeros(numrows,(m+1)^n);
dimensionmatrix = zeros(numrows,1);
bestdimensionmatrix = zeros(numrows,1);
dimdiffmatrix = zeros(numrows,1);
classicalboundmatrix = zeros(numrows,1);
quantumboundmatrix = zeros(numrows,1);
beststatematrix = zeros(numrows,8);

fprintf("Starting...\n")

for row = 1:numrows
    ketstate = states(row,:).';
    b = rhodecomp(ketstate * ketstate').';
    numelx = size(tensorprojs,3);
    A = vectensorprojs;

    t1 = clock;
    
    cvx_begin quiet; 
        variable x(numelx);
        dual variable y;
        minimize( norm( x, 1 ) );
        subject to;
          y: A * x == b;
    cvx_end;
    robustness=cvx_optval;
    
    [dimension,smax] = calcdimandclassicalbound(n,d,m,y.');
    
    robustnessmatrix(row,:) = norm(x,1);
    inequalitymatrix(row,:) = y.';
    dimensionmatrix(row,:) = dimension;
    bestdimensionmatrix(row,:) = bestdimension;
    dimdiffmatrix(row,:) = dimension - bestdimension;
    classicalboundmatrix(row,:) = smax;
    
    B = 0;
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
        B = B + inequalitymatrix(row,index)*kron(tensprod,oparray(:,:,3));
    end
    
    [matrixeigvec,matrixeig] = eig(B);
    
    eigenvalues = zeros(1,length(matrixeig));
    for i1 = 1:length(matrixeig)
    eigenvalues(i1) = matrixeig(i1,i1);
    end
    quantum = max(eigenvalues);
    index = find(eigenvalues == quantum);
    beststate = matrixeigvec(:,index);
   
    quantumboundmatrix(row,:) = quantum;
    beststatematrix(row,:) = beststate.';       
    
    T.Properties.RowNames = row;
    
    % Get the current time in seconds.
    t2 = clock;
    % Calculate the time taken for this calculation.
    numsecs = etime(t2,t1);
    % Add this to an array.
    numsecsarray(row) = numsecs;
    % Estimated time left in seconds, from the average of the time taken per
    % calculation.
    timeleftsecs = (numrows-row)*mean(numsecsarray);
    % Print out the estimated time left in a nice format.
    estimatedtimeleft = datestr(timeleftsecs/(24*60*60), 'DD:HH:MM:SS.FFF')  
    
end

data = table(states,robustnessmatrix,inequalitymatrix,dimensionmatrix,bestdimensionmatrix,dimdiffmatrix,classicalboundmatrix,quantumboundmatrix,beststatematrix);
