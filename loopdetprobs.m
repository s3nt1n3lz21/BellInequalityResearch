function [smax,varsgivesmax,amountofvarsgivesmax] = loopdetprobs(varstoloop,vars,coefflist,n,d,m)
% LOOPDETPROBS Loop over all the possible values of the single party deterministic probability distributions (deterministic strategies) and calculate the value of the bell inequality through recursion

% smax is the current branch's maximum value of the expression so far
% varsgivesmax is a matrix whoose rows are the values of the deterministic probabilities (0 or 1) that give the current branch's smax
% amountofvarsgivesmax is the number of rows in varsgivesmax, this is (needed)?

    varsgivesmax = NaN(2^(n*m*d),n*d*m);
    smax = 'x';
    amountofvarsgivesmax = 0;
    if varstoloop >= 1
    % if there is still a variable to loop over, loop over it's possible values (0 or 1)            
         for value = 0:1
            vars(varstoloop) = value;
            % then loop over the rest by calling the function again, giving it the current values of the variables (var)          
            [innersmax,innervarsgivesmax,inneramountofvarsgivesmax] = loopdetprobs(varstoloop-1,vars,coefflist,n,d,m);
            if smax == 'x'
                if innersmax == 'x'
                    % do nothing
                else
                    % if smax still undefined but innersmax has a value set the value of smax                    
                    smax = innersmax;
                    varsgivesmax(1:inneramountofvarsgivesmax,:) = innervarsgivesmax(1:inneramountofvarsgivesmax,:);
                    amountofvarsgivesmax = inneramountofvarsgivesmax;
                end    
            else
                if innersmax == 'x'
                    % do nothing                     
                else
                    if innersmax > smax
                        smax = innersmax;
                        % reset the matrix of vars that gives smax with these values                
                        varsgivesmax = NaN(2^(n*m*d),n*m*d);
                        innervarsgivesmax(1:inneramountofvarsgivesmax,:);
                        varsgivesmax(1:inneramountofvarsgivesmax,:) = innervarsgivesmax(1:inneramountofvarsgivesmax,:);
                        amountofvarsgivesmax = inneramountofvarsgivesmax;
                    elseif innersmax == smax
                        % append the values                    
                        varsgivesmax(amountofvarsgivesmax+1:amountofvarsgivesmax+inneramountofvarsgivesmax,:) = innervarsgivesmax(1:inneramountofvarsgivesmax,:);
                        amountofvarsgivesmax = amountofvarsgivesmax + inneramountofvarsgivesmax;
                    else
                        % do nothing                     
                    end
                end
            end
         end
    else
    % If no more variables to loop over calculate the expression for the current values of the variables
        % check the normalization/no signalling contraints to the values of these variables
        % If the current values of the variables obey the conditions then continue, else skip the calculation
        % Use a flag "constraintsobeyed" to check whether they are obeyed, initialise this to true        
        constraintsobeyed = true;
        % for each party ni and setting mi check the sum of the outcomes is unity. i.e Sum over d_i of D_ni(d_i|m_i) = 1 
        for ni = 1:n
            % if at any point the constaints are not obeyed then stop carrying on to check they are obeyed           
            if not(constraintsobeyed)
                break;
            else
                for mi = 1:m
                    % if at any point the constaints are not obeyed then stop carrying on to check they are obeyed 
                    if not(constraintsobeyed)
                        break;
                    else
                        % check that the constaints are obeyed                        
                        constraintexpression = 0;
                        for di = 1:d
                            constraintexpression = constraintexpression + vars(f(ni,di,mi,d,m));
                        end
                        if not(constraintexpression == 1)
                            constraintsobeyed = false;
                        end
                    end
                end
            end
        end
        if constraintsobeyed
            smax = 0;
            % for each term in the corrcoefflist determine which correlators we have to calculate            
            for i1 = 1:numel(coefflist)
                coeff = coefflist(i1);
                % if the coefficient is zero skip the calculation 
                if not(coeff) == 0
                    % for each correlator e.g <m1 m2> calculate the contributing terms to the expression for s              
                    % get the correct settings e.g m1 and m2 from how far into the vector coefflist we are           
                    marray = getsettings(i1,m,n);
                    % find which of the parties do make measurements as this affects the form of the expression                    
                    partiesmakemeasurements = find(marray);
                    % from this find how many variables we have to loop/sum over                    
                    numdvars = length(partiesmakemeasurements);
                    % create an array to hold these values as they change                    
                    darray = zeros(1,n);
                    % now calculate the correlator of measurements e.g m1 and m2 by looping with recursion                    
                    finalcurrentcorr = calcinnercorr(numdvars,marray,darray,vars,d,m,partiesmakemeasurements);                
                    smax = smax + coeff*finalcurrentcorr;
                end
            end
            varsgivesmax(1,:) = vars;
            amountofvarsgivesmax = 1;
        else
        % if constaints not obeyed make s a character 'x' to say there isn't a value for s yet (NaN always returns false in conditional statements)       
            smax = 'x';
            amountofvarsgivesmax = 0;
            varsgivesmax = NaN(2^(n*m*d),n*m*d);
        end
    end