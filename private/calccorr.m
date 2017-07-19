function currcorr = calccorr(varstoloop,marray,darray,detprobs,d,m,partiesmakemeasurements)
% function to calculate the correlator of a set of measurements m1 m2 ..
    % if there are still variables (outcomes) to loop over then loop over them    
    if varstoloop >= 1
        corsofar = 0;
        % for di = 1:d         
        for dvalue = 1:d    
            % only loop over the outcomes of the parties that take measurements             
            darray(partiesmakemeasurements(varstoloop)) = dvalue;
            innercorr = calccorr(varstoloop-1,marray,darray,detprobs,d,m,partiesmakemeasurements);
            corsofar = corsofar + innercorr;
        end
        currcorr = corsofar;
    % if there are no more to loop over then calculate the contribution to the correlator    
    else
        curprodterm = 1;
        % the expression will take the form of k products if there are k parties that make measurements        
        for party = partiesmakemeasurements
            curprodterm = curprodterm*detprobs(getdetprobindex(party,darray(party),marray(party),d,m));
        end
        % also multiply by the prefactor        
        curprodterm = curprodterm*((-1)^(sum(darray)));
        % return this value to the parent loop or final expression        
        currcorr = curprodterm;
    end