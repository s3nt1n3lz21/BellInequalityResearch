% Main file

paulis = cat(3,[1,0;0,1],[0,1;1,0],[0,-i;i,0],[1,0;0,-1]);

for j = 1:3
    sproj(:,:,j) = 0.5*(eye(2)+paulis(:,:,j+1));
    sproj(:,:,j+3) = 0.5*(eye(2)-paulis(:,:,j+1));
end

sproj

ind = 1;
for j = 1:6
    for k = 1:6
        PProj(:,:,ind) = kron(sproj(:,:,j),sproj(:,:,k)) 
        vecPProj(:,ind) = rhodecomp(kron(sproj(:,:,j),sproj(:,:,k)))
        ind = ind+1;
    end
end

PProj
vecPProj

states = importdata('Three_qubit_hierarchy.mat');
ketstate = states(1,:).'
rho = ketstate*ketstate'
b = rhodecomp(rho)

% for row = 1:size(states,1)
%     ketstate = states(row,:).'
%     b=rhodecomp(ketstate * ketstate').';
%     n=size(PProj,3);
%     A=vecPProj;
%     
%     cvx_begin 
%         variable x(n)
%         dual variable y;
%         minimize( norm( x, 1 ) )
%         subject to
%           y: A * x == b
%     cvx_end
%     robustness=cvx_optval;
%     
%     y
% end