classdef cbanddimcalc < handle
% CBANDDIMCALC A class used to calculate the dimension and classical bound
% of a bell inequality.
   properties      
      % The number of parties n, measurement outcomes d and measurement settings m. Here they are assumed to all be the same.      
      n
      d
      m
      % The number of single party deterministic probabilities.
      numvars
      % The maximum possible dimension of the bell inequality is an estimate of
      % the maximum number of ways of getting smax, and so the memory allocated.
      maxdim
      % The correlator coefficient list.  
      corrcoefflist
      % The array of the current values of the single party deterministic probabilities.      
      detprobvalues
      % The current maximum found so far.
      smax
      % The array of deterministic probabilities that give smax.      
      detprobsgivesmax
      % How many rows full is the array of deterministic probabilities (not the number of rows).     
      detprobsrows
      % An array to hold the current measurement settings.      
      mvalues
      % An array used in calccorr to hold the values of the measurement outcomes d for each party as it loops over them.     
      dvalues
      % An array used in calccorr to hold all of the correlator terms which it
      % calculates.
      corrvalues
      % current correlator term for the (Parties that Make Measurements).
      pmm
      % How full is the array of corrvalues.
      corrvaluesrows
      % An array used in calcprobdists to hold the current values of the possible
      % outcomes.
      darray
      % An array used in calcprobdists to hold the current values of the
      % measurent settings.
      marray
      % An array used in calcprobdists to hold the current probability
      % distribution.
      detprobs
      % An array used to hold the current probability distribution in
      % calcprobdists
      probdist
      % How full is the vector probdist
      probdistelements
   end
   methods       
      function obj = cbanddimcalc(n,d,m,corrcoefflist)
         % CBANDDIMCALC Constructor function to the initialise properties.      
         obj.maxdim = ((m*(d-1)+1)^n) - 1;
         obj.n = n;
         obj.d = d;
         obj.m = m;
         obj.numvars = n*m*d;
         obj.corrcoefflist = corrcoefflist;
         obj.detprobvalues = zeros(1,obj.numvars);
         obj.detprobsgivesmax = zeros(obj.maxdim,n*m*d);
         obj.detprobsrows = 0;
         obj.smax = 'x';
         obj.dvalues = zeros(1,n);
      end
      function [dim,smax] = calc(obj)
      % CALC Calculate the dimension and classical bound.
       loopdetprobs(obj,obj.numvars);
       % Calculate the probability distribution vectors from determinsitic
       % probability array.
       probdistsgivesmax = calcprobdists(obj);
       % Calculate the dimension.
       dim = calcdim(obj,probdistsgivesmax);
       smax = obj.smax;
      end
      function loopdetprobs(obj,varstoloop)
      % LOOPDETPROBS The main algorithm to loop over the possible values of the
      % deterministic probabilities (0 or 1).
      
          % If there is a still a deterministic probability to loop over then loop
          % over its possible values.
          if varstoloop >= 1
               for value = 0:1
                  obj.detprobvalues(varstoloop) = value;
                  % Then loop over the rest by calling the function again, with one less variable to loop over.    
                  loopdetprobs(obj,varstoloop-1);
               end
               
          % If there are no more deterministic probabilities to loop over then
          % calculate the Bell Inequality correlator for the current values of the
          % deterministic probabilities.
          else

              % First check that for this particular arrangement of probabilities, the normalization and no-signalling 
              % contraints are obeyed. Use a flag "constraintsobeyed" to check whether they are obeyed, initialise this to true.        
              constraintsobeyed = true;
              % For each possible choice of ni and mi check the sum of the outcomes is unity. i.e The sum over di of D_ni(di|mi) = 1 
              for ni = 1: obj.n
                  % If at any point the constaints are not obeyed then stop carrying on to check that they are obeyed.           
                  if not(constraintsobeyed)
                      break;
                  else
                      for mi = 1: obj.m
                          % If at any point the constaints are not obeyed then stop carrying on to check that they are obeyed.
                          if not(constraintsobeyed)
                              break;
                          else
                              % Check that the constaints are obeyed.                        
                              constraintexpression = 0;
                              for di = 1: obj.d
                                  constraintexpression = constraintexpression + obj.detprobvalues(getdetprobindex(obj,ni,di,mi));
                              end
                              if not(constraintexpression == 1)
                                  constraintsobeyed = false;
                              end
                          end
                      end
                  end
              end
            % If the constraints are obeyed then continue with the calculation.            
            if constraintsobeyed
                % Initialise the value of the expression to be 0.                
                s = 0;
                % For each term in the correlator coefficient list calculate the corresponding correlator. e.g <m1 m2>          
                for i1 = 1:numel(obj.corrcoefflist)
                    coeff = obj.corrcoefflist(i1);
                    % If the coefficient is zero then just skip the calculation otherwise continue. 
                    if not(coeff) == 0           
                        % Get the correct measurement settings (e.g m1 and m2) from how far into
                        % the correlator coefficient list we are.
                        obj.mvalues = dec2base(i1-1,obj.m+1,obj.n)-'0';
                        % Find which of the parties do make measurements as this affects the form of the expression to be calculated. 
                        % If k parties make measurements then there will be k products of variables.            
                        obj.pmm = find(obj.mvalues);
                        % The NUMber of Parties that Make Measurements is the number of outcome variables di we have to loop over.              
                        numpmm = length(obj.pmm);

                        % Initialise the outcome variables, the array of correlator values and numbers of rows it is full before each calculation.  
                        obj.dvalues = zeros(1,obj.n);
                        obj.corrvalues = zeros((obj.d)^(numpmm),1);
                        obj.corrvaluesrows = 0;
                        
                        % Now calculate the correlator of the measurements (e.g m1 and m2) by looping with recursion.                       
                        calccorr(obj,numpmm);
                        
                        % Once all of the correlator values have been calculated calculate the sum and then add this to value of the expression 
                        % first multiplying by the corresponding coefficient in the correlator coefficient list.                        
                        currentcorr = sum(obj.corrvalues);
                        s = s + coeff*currentcorr;
                    end
                end
                % If smax is not yet defined then set the value of smax and the array of deterministic probabilities that give smax. 
                if obj.smax == 'x'
                  obj.smax = s;
                  obj.detprobsgivesmax(1,:) = obj.detprobvalues;
                  obj.detprobsrows = obj.detprobsrows + 1;
                else
                     % If s is greater than the current maximum then set the new value of smax and reset the array of deterministic 
                     % probabilities that give smax with this new value.
                     if s > obj.smax
                       obj.smax = s;
                       obj.detprobsgivesmax = zeros(obj.maxdim,obj.numvars);
                       obj.detprobsgivesmax(1,:) = obj.detprobvalues;
                       obj.detprobsrows = 1;

                     % If s is the same as the current maximum, append the deterministic
                     % probabilities to the array of deterministic probabilities that give smax.
                     elseif s == obj.smax
                       obj.detprobsgivesmax(obj.detprobsrows+1,:) = obj.detprobvalues;
                       obj.detprobsrows = obj.detprobsrows + 1;
                       
                     % Else If s is less than smax don't do anything, just continue to the next loop.          
                     else
                      
                     end    
                end
                
            % If the constaints are not obeyed for this particular combination of deterministic probabilities 
            % then don't do anything, just continue to the next loop.
            else

            end
          end
      end
      function index = getdetprobindex(obj,ni,di,mi)
          % GETDETPROBINDEX Calculates the index of the deterministic probability with D_ni(di,mi) within the array of deterministic probabilities.          
          index = (ni-1)*obj.d*obj.m+(di-1)*obj.m+mi; 
      end
      function calccorr(obj,varstoloop)
      % CALCCORR Calculates the correlator of a set of measurements (e.g <m1 m2>).
      
      % If there are still outcome variables to loop over, then loop over them. Only loop over the possible outcomes
      % of the parties that make measurements.
       if varstoloop >= 1   
           for dvalue = 1: obj.d              
               obj.dvalues(obj.pmm(varstoloop)) = dvalue;
               % Then call the function again with one less variable to loop over.               
               calccorr(obj,varstoloop-1);
           end
       % If there are no more to loop over then calculate the current contribution to this correlator.    
       else
           curprodterm = 1;
           % The expression will take the form of k products of probabilities if there
           % are k parties that make measurements. For each party that does make
           % measurements calculate the corresponding product.
           for party = obj.pmm
               curprodterm = curprodterm*obj.detprobvalues(getdetprobindex(obj,party,obj.dvalues(party),obj.mvalues(party)));
           end
           % Then multiply by the prefactor. Note the other dvalues will stay zero and so not contribute.        
           curprodterm = curprodterm*((-1)^(sum(obj.dvalues)));
           % Add this contribution to the correlator to the array corrvalues and keep track how full the array is.  
           obj.corrvalues(obj.corrvaluesrows+1) = curprodterm;
           obj.corrvaluesrows = obj.corrvaluesrows + 1;
       end
      end
      function [probdistsgivesmax] = calcprobdists(obj)
      % CALCPROBDISTS Calculate the probability distribution vectors that give
      % the classical bound from the array of deterministic probabilities that give the classical bound.
          obj.marray = zeros(1,obj.n);
          obj.darray = zeros(1,obj.n);
          probdistsgivesmax = zeros(obj.detprobsrows,(obj.d*obj.m)^obj.n);
          for row = 1:obj.detprobsrows
              obj.detprobs = obj.detprobsgivesmax(row,:);
              obj.probdist = NaN(1,(obj.d*obj.m)^obj.n);
              obj.probdistelements = 0;
              calcprobdist(obj,obj.n,obj.n);
              probdistsgivesmax(row,:) = obj.probdist;
          end
      end
      function calcprobdist(obj,dvarstoloop,mvarstoloop)
      % CALCPROBDIST Calculates the probability distribution vector from the
      % deterministic probabilities.
      
          % If there are outcome variables to loop over then loop over them (e.g d1
          % d2).
          if dvarstoloop >= 1
              for dvalue = 1:obj.d
                  obj.darray(dvarstoloop) = dvalue;
                  calcprobdist(obj,dvarstoloop-1,mvarstoloop);
              end
          % If there are measurement settings to loop over then loop over them (e.g
          % m1 m2).
          elseif mvarstoloop >= 1
              for mvalue = 1:obj.m
                  obj.marray(mvarstoloop) = mvalue;
                  calcprobdist(obj,dvarstoloop,mvarstoloop-1);
              end
          % Otherwise calculate the current probability term in the probability distribution vector.          
          else
              curprodterm = 1;
              for party = 1:obj.n
                  curprodterm = curprodterm*obj.detprobs(getdetprobindex(obj,party,obj.darray(party),obj.marray(party)));
              end
              obj.probdist(obj.probdistelements+1) = curprodterm;
              obj.probdistelements = obj.probdistelements + 1;
          end
      end
      function [dim] = calcdim(~,probdistsgivesmax)
      % CALCDIM Calculate the dimension of the bell inequality from the
      % probability distributions that give the classical bound.
      reducedmatrix = rref(probdistsgivesmax);
      dim = rank(reducedmatrix);
      end
   end
end