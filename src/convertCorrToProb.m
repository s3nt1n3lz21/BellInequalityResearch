function [maxNoMeasOutcomesList,probCoeffList] = convertCorrToProb(nvar,d,m,corrcoefflist)

    for i1 = 1:length(corrcoefflist)
        msettings = str2num(dec2base(i1-1,m+1,nvar)')';
        numpmm = nnz(msettings);
        pmm = find(msettings);
        if corrcoefflist(i1)
            currentvec = zeros(1,d^(numpmm));
            for i2 = 1:d^(numpmm)
                dvalues = str2num(dec2base(i2-1,d,numpmm)')'+ones(1,numpmm);
                currentvec(i2) = corrcoefflist(i1)*(-1)^(sum(dvalues)-numpmm);
            end
        else
            currentvec = zeros(1,d^(numpmm));
        end
        probCoeffList(:,i1) = {currentvec'};
    end
        
    dlistelement = (ones(1,m)*d);
    for i3 = 1:nvar
        maxNoMeasOutcomesList(:,i3) = {dlistelement'};
    end