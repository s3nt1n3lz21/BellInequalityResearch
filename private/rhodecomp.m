function rhodecomp=rhodecomp(rho)

n=log2(length(rho));
paulis=cat(3,[1,0;0,1],[0,1;1,0],[0,-i;i,0],[1,0;0,-1])
for j=0:4^n-1
    decomp(j+1)=real(trace(kronmult(mat2cell(paulis(:,:,dec2basevec(j,4,n)+1), 2, 2, ones(1,n)),rho)));
end

rhodecomp=decomp;