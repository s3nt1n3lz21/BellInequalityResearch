function [probdistsofar,numelements] = calcprobdist(vars,dvarstoloop,mvarstoloop,darray,marray,n,d,m)
% CALCPROBDIST Calculates the probability distribution vector from the deterministic probabilities
    probdistsofar = NaN(1,(d*m)^n);
    if dvarstoloop >= 1
        numelements = 0;
        for dvalue = 1:d
            darray(dvarstoloop) = dvalue;
            [innerelements,numinnerelements] = calcprobdist(vars,dvarstoloop-1,mvarstoloop,darray,marray,n,d,m);
            probdistsofar(numelements+1:numelements+numinnerelements) = innerelements(1:numinnerelements);
            numelements = numelements+numinnerelements;
        end
    elseif mvarstoloop >= 1
        numelements = 0;
        for mvalue = 1:m
            marray(mvarstoloop) = mvalue;
            [innerelements,numinnerelements] = calcprobdist(vars,dvarstoloop,mvarstoloop-1,darray,marray,n,d,m);
            probdistsofar(numelements+1:numelements+numinnerelements) = innerelements(1:numinnerelements);
            numelements = numelements+numinnerelements;
        end
    else
        numelements = 1;
        curprodterm = 1;
        for party = 1:n
            curprodterm = curprodterm*vars(getdetprobindex(party,darray(party),marray(party),d,m));
        end
        probdistsofar(numelements) = curprodterm;
    end