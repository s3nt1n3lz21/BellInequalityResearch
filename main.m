% Here I calculate the robustness and other quantities for a variety of
% states in the file "three_qubit_hierarchary". I calculate the classical
% bound and dimension of the bell inequality found by calculating the
% robustness. I also calculate the difference of the dimension of the bell
% inequality from the dimension of the best possible bell inequality and an
% estimation of the quantum bound using the eigenvalue of the bell
% operator. The corresponding eigenvector represents a better quantum state
% that can be used to test the bell inequality.

% define the pauli operators
paulis = cat(3,[1,0;0,1],[0,1;1,0],[0,-i;i,0],[1,0;0,-1]);

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

states = importdata('Three_qubit_hierarchy.mat');
numrows = size(states,1);
ymatrix = zeros(numrows,(m+1)^n);
for row = 1:1
    ketstate = states(row,:).';
    b = rhodecomp(ketstate * ketstate').';
    numelx = size(tensorprojs,3);
    A = vectensorprojs;
    t1 = clock;
    
    cvx_begin 
        variable x(numelx)
        dual variable y;
        minimize( norm( x, 1 ) )
        subject to
          y: A * x == b
    cvx_end
    robustness=cvx_optval;
    
    ymatrix(row,:) = y.';
    y.'
    [dimension,smax] = calcdimandclassicalbound(n,d,m,ymatrix(row,:))
    
    t2 = clock;
    numsecs = etime(t2,t1)
    timeleft = (numrows-row)*numsecs
    
end

