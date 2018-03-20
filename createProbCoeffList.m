function [probCoeffList] = createProbCoeffList(nvar,d,m,ylist)

    yLeftOver = ylist;
    for i1 = 1:(m+1)^nvar
        msettings = str2num(dec2base(i1-1,m+1,nvar)')';
        numpmm = nnz(msettings);
        if numpmm == length(msettings)
            ySlice = yLeftOver(1:2^d);
            yLeftOver = yLeftOver(2^d+1:length(yLeftOver));
            currentvec = ySlice;
        else    
            currentvec = zeros(1,d^(numpmm));
        end
        probCoeffList(:,i1) = {currentvec'};
    end