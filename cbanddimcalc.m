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
      % calculates
      corrvalues
      % current correlator term for the (Parties that Make Measurements)
      pmm
      % How full is the array of corrvalues
      corrvaluesrows
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
       dim = obj.dim;
       smax = obj.smax;
      end
      function loopdetprobs(obj,varstoloop)
      % LOOPDETPROBS The main algorithm to loop over the possible values of the deterministic probabilities (0 or 1)       
          if varstoloop >= 1
          % if there is still a variable to loop over, loop over it's possible values (0 or 1)            
               for value = 0:1
                  obj.detprobvalues(varstoloop) = value;
                  % then loop over the rest by calling the function again    
                  loopdetprobs(obj,varstoloop-1);
               end
          else
          % If no more variables to loop over calculate the expression for the current values of the variables
              % check the normalization/no signalling contraints to the values of these variables
              % If the current values of the variables obey the conditions then continue, else skip the calculation
              % Use a flag "constraintsobeyed" to check whether they are obeyed, initialise this to true        
              constraintsobeyed = true;
              % for each party ni and setting mi check the sum of the outcomes is unity. i.e Sum over d_i of D_ni(d_i|m_i) = 1 
              for ni = 1: obj.n
                  % if at any point the constaints are not obeyed then stop carrying on to check they are obeyed           
                  if not(constraintsobeyed)
                      break;
                  else
                      for mi = 1: obj.m
                          % if at any point the constaints are not obeyed then stop carrying on to check they are obeyed 
                          if not(constraintsobeyed)
                              break;
                          else
                              % check that the constaints are obeyed                        
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
            if constraintsobeyed
                currentsmax = 0;
                % for each term in the corrcoefflist determine which correlators we have to calculate            
                for i1 = 1:numel(obj.corrcoefflist)
                    coeff = obj.corrcoefflist(i1);
                    % if the coefficient is zero skip the calculation 
                    if not(coeff) == 0
                        % for each correlator e.g <m1 m2> calculate the contributing terms to the expression for s              
                        % get the correct settings e.g m1 and m2 from how far into the vector coefflist we are
                        obj.mvalues = dec2base(i1-1,obj.m+1,obj.n)-'0';
                        % find which of the parties do make measurements as this affects the form of the expression                    
                        obj.pmm = find(obj.mvalues);
                        % from this find how many variables we have to loop/sum over                    
                        numpmm = length(obj.pmm);

                        % Reset the dvalues, corrvalues and rows before each calculation    
                        obj.corrvaluesrows = 0;
                        obj.corrvalues = zeros((obj.d)^(numpmm),1);
                        obj.dvalues = zeros(1,obj.n);
                        
                        % now calculate the correlator of measurements e.g m1 and m2 by looping with recursion                       
                        calccorr(obj,numpmm);
                        currentcorr = sum(obj.corrvalues);
                        currentsmax = currentsmax + coeff*currentcorr;
                    end
                end
                if obj.smax == 'x'
                  % if smax not yet defined then set the value of smax and array of deterministic probabilities that give smax                    
                  obj.smax = currentsmax;
                  obj.detprobsgivesmax(1,:) = obj.detprobvalues;
                  obj.detprobsrows = obj.detprobsrows + 1;
                else
                     if currentsmax > obj.smax
                     % Set the new value of smax and reset the array of deterministic probabilities that give smax with this new value
                       obj.smax = currentsmax;
                       obj.detprobsgivesmax = zeros(obj.maxdim,obj.numvars);
                       obj.detprobsgivesmax(1,:) = obj.detprobvalues;
                       obj.detprobsrows = 1;
                     elseif currentsmax == obj.smax
                     % Append the values of detprobs
                       obj.detprobsgivesmax(obj.detprobsrows+1,:) = obj.detprobvalues;
                       obj.detprobsrows = obj.detprobsrows + 1;
                     else
                     % otherwise don't do anything, continue to the next loop of deterministic probabilities                    
                     end    
                end
            else
            % if constaints not obeyed don't do anything, just continue to the next loop of deterministic probabilities    
            end
          end
      end
      function index = getdetprobindex(obj,ni,di,mi)
        index = (ni-1)*obj.d*obj.m+(di-1)*obj.m+mi; 
      end
      function calccorr(obj,varstoloop)
% function to calculate the correlator of a set of measurements m1 m2 ..
    % if there are still variables (outcomes) to loop over then loop over them    
       if varstoloop >= 1
           % for di = 1:d         
           for dvalue = 1: obj.d    
               % only loop over the outcomes of the parties that take measurements             
               obj.dvalues(obj.pmm(varstoloop)) = dvalue;
               calccorr(obj,varstoloop-1);
           end
       % if there are no more to loop over then calculate the contribution to the correlator    
       else
           curprodterm = 1;
           % the expression will take the form of k products if there are k parties that make measurements
           for party = obj.pmm
               curprodterm = curprodterm*obj.detprobvalues(getdetprobindex(obj,party,obj.dvalues(party),obj.mvalues(party)));
           end
           % also multiply by the prefactor        
           curprodterm = curprodterm*((-1)^(sum(obj.dvalues)));
           % Add this value to the array corrvalues and keep track how full the array is  
           obj.corrvalues(obj.corrvaluesrows+1) = curprodterm;
           obj.corrvaluesrows = obj.corrvaluesrows + 1;
       end
      end
      function calcprobdists(obj)
          darray = zeros(1,obj.n);
          marray = zeros(1,obj.n);
          probdistsgivesmax = zeros(numdetprobs,(obj.d*obj.m)^obj.n);
          for row = 1:obj.detprobsrows
              detprobs = obj.detprobsgivesmax(row,:);
              probdist = calcprobdist(detprobs,obj.n,obj.n,darray,marray,n,d,m);
              probdistsgivesmax(row,:) = probdist;
          end
          reducedmatrix = rref(probdistsgivesmax)
          dimension = rank(reducedmatrix);
      end
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
      end
   end
end