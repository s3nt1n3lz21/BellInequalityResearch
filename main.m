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
tensorprojs = zeros(8,8,64);
vectensorprojs(:
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

for row = 1:size(states,1)
    ketstate = states(row,:).'
    b=rhodecomp(ketstate * ketstate').'
    lb = length(b)
    n=size(tensorprojs,3);
    A=vectensorprojs
    sizeA = size(A)
    
    cvx_begin 
        variable x(n)
        dual variable y;
        minimize( norm( x, 1 ) )
        subject to
          y: A * x == b
    cvx_end
    robustness=cvx_optval;
    
    y
end