classdef classicalboundcalculator
   properties
      % The number of parties, measurement outcomes and measurement settings       
      n
      d
      m
      % The number of deterministic probability distributions
      numvars
      % The maximum possible dimension of the bell inequality is also the maximum number of deterministic probability distributions that give smax
      maxdim
      % The correlator coefficient list  
      corrcoefflist
      % The array of current values of the deterministic probabilities      
      detprobvalues
      % The current maximum found so far
      smax
      % The array of deterministic probabilities that give smax      
      detprobsgivesmax
      % How full is the array of deterministic probabilities      
      detprobsrows
      % An array to hold the current measurement settings      
      mvalues
      % An array used in calccorr to loop over the possible values of the measurement outcomes      
      dvalues
      % An array used in calccorr to hold all of the correlator terms which it calculates and then is summed later
      corrvalues
      % An array used in calccorr to hold the measurement settings for the
      % current correlator term for the (Parties that Make Measurements)
      pmm
      % How full is the array of corrvalues
      corrvaluesrows
   end
   methods       
      function obj = classicalboundcalculator(n,d,m,corrcoefflist)
         % Constructor function to initialise properties      
         obj.maxdim = 2*(d-1)*m+((d-1)^2)*m^2;
         obj.n = n;
         obj.d = d;
         obj.m = m;
         obj.numvars = n*m*d;
         obj.corrcoefflist = corrcoefflist;
         obj.detprobvalues = zeros(obj.numvars);
         obj.detprobsgivesmax = zeros(obj.maxdim,n*m*d);
         obj.detprobsrows = 0;
         obj.smax = 'x';
         obj.dvalues = zeros(1,n);
      end
      function [smax,detprobsgivesmax] = calcclassicalbound(obj)
       loopdetprobs(obj,obj.numvars);
       smax = obj.smax;
       detprobsgivesmax = obj.detprobsgivesmax;
      end
      function loopdetprobs(obj,varstoloop)
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
                        obj.mvalues = getsettings(obj,i1);
                        % find which of the parties do make measurements as this affects the form of the expression                    
                        obj.pmm = find(obj.mvalues);
                        % from this find how many variables we have to loop/sum over                    
                        numpmm = length(obj.pmm);
                        % now calculate the correlator of measurements e.g m1 and m2 by looping with recursion
                        obj.corrvaluesrows = 0;
                        obj.corrvalues = zeros((obj.d)^(numpmm),1);
                        calccorr(obj,numpmm);
                        currentcorr = sum(obj.corrvalues);
                        currentsmax = currentsmax + coeff*currentcorr;
                    end
                end
                if obj.detprobsrows+1 > size(obj.detprobsgivesmax,1)
                  % If the index has gone out of bounds then we know that the current value of smax cannot be the correct one so reset the array
                   obj.detprobsgivesmax = zeros(obj.maxdim,obj.numvars);
                   obj.detprobsrows = 0;
                end
                if obj.smax == 'x'
                  % if smax not yet defined then set the value of smax and array of deterministic probabilities that give smax                    
                  obj.smax = currentsmax;
                  obj.detprobsgivesmax(obj.detprobsrows+1,:) = obj.detprobvalues;
                  obj.detprobsrows = obj.detprobsrows + 1;
                else
                  if currentsmax > obj.smax
                  % Set the new value of smax and reset the array of deterministic probabilities that give smax                    
                    obj.smax = currentsmax;
                    obj.detprobsgivesmax = zeros(obj.maxdim,obj.numvars);
                    obj.detprobsrows = 0;
                  elseif currentsmax == obj.smax
                  % Append the values of detprobs
                    obj.detprobsgivesmax(obj.detprobsrows+1,:) = obj.detprobvalues;
                    obj.detprobsrows = obj.detprobsrows + 1;
                  else
                  % otherwise don't do anything, continue to the next loop                     
                  end    
                end
            else
            % if constaints not obeyed don't do anything, just continue to the next loop     
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
               curprodterm = curprodterm*obj.detprobs(obj.getdetprobindex(obj,party,darray(party),marray(party)));
           end
           % also multiply by the prefactor        
           curprodterm = curprodterm*((-1)^(sum(obj.dvalues)));
           % Add this value to the array corrvalues and keep track how full the array is  
           obj.corrvalues(obj.corrvaluesrow+1) = curprodterm;
           obj.corrvaluesrow = obj.corrvaluesrow + 1;
       end
      end
      function settings = getsettings(obj,index)
       % function to calculate what the measurements settings are for the current correlator
       % create an array to hold these values    
       m = obj.m;
       n = obj.n;
       settings = zeros(1,n);
       % in general the index of a particular group of settings is given by "index = 1 + m1*(m+1)^(n-1) + m2*(m+1)^(n-2) + m3*(m+1)^(n-3)+..."    
       settings(1) = idivide(int32(index-1),(m+1)^(n-1));
       remainder = (index-1)/((m+1)^(n-1))-settings(1);
       for i1 = 2:n
           settings(i1) = round(remainder*(m+1));
           remainder = remainder*(m+1) - settings(i1);
       end
      end
   end
end