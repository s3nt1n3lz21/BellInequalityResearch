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
      % The array to hold the di values of each constraint equation that are 1.
      indexarray
      % The maximum possible dimension of the bell inequality is an estimate of
      % the maximum number of ways of getting smax, and so the memory allocated.
      maxdim
      % The correlator coefficient list.  
      coefflist
      % The array of the current values of the single party deterministic probabilities.      
      detprobvalues
      % The current maximum found so far.
      smax
      % The array of deterministic probabilities that give smax.      
      detprobsgivesmax
      % How many rows full is the array of deterministic probabilities (not the number of rows).     
      detprobsrows
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
      function obj = cbanddimcalc(n,d,m,coefflist)
         % CBANDDIMCALC Constructor function to the initialise properties.      
         obj.maxdim = ((m*(d-1)+1)^n) - 1;
         obj.n = n;
         obj.d = d;
         obj.m = m;
         obj.numvars = n*m*d;
         obj.coefflist = coefflist;
         obj.detprobvalues = zeros(1,obj.numvars);
         %Make detprobsgivesmax dynamic but specify an upper bound.
         obj.detprobsgivesmax = zeros(obj.d^(obj.n*obj.m),obj.numvars);
         coder.varsize('obj.detprobsgivesmax',[obj.d^(obj.n*obj.m),obj.numvars]);
         obj.detprobsrows = 0;
         obj.smax = 'x';
         obj.indexarray = zeros(1,(n*m));
         %listsize = size(coefflist,1);
         %if not(listsize == (m+1)^n)
         %   fprintf("Error, the dimension of the correlator coefficient list does not match the scenario\n The dimension is %.0f when it should be %.0f",size(corrcoefflist,2),(m+1)^n)
         %end
      end
      function [dim,smax] = calc(obj)
      % CALC Calculate the dimension and classical bound.
       loopdetprobs(obj,obj.n*obj.m);
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
               for value = 1:obj.d
                  obj.indexarray(varstoloop) = value;
                  % Then loop over the rest by calling the function again, with one less variable to loop over.    
                  loopdetprobs(obj,varstoloop-1);
               end
               
          % If there are no more deterministic probabilities to loop over then
          % calculate the Bell Inequality correlator for the current values of the
          % deterministic probabilities.
          else

                % Initialise the value of the expression to be 0 and get the values of the variables from the indexes.         
                obj.detprobvalues = zeros(1,obj.numvars);
                for ni = 1:obj.n
                    for mi = 1:obj.m
                    dindex =  (ni-1)*obj.m+mi;
                    di = obj.indexarray(dindex);
                    obj.detprobvalues(getdetprobindex(obj,ni,di,mi)) = 1;
                    end
                end

                s = 0;
                % For each term in the coefficient list calculate the corresponding contribution from that probability. e.g coeff*prob    
                for i2 = 1:size(obj.coefflist,1)
                    % Get the correct measurement settings (e.g m1 and m2) from how far into the coefficient list we are.
                    mvalues = dec2base(i2-1,obj.m+1,obj.n)-'0';
                    % Find which of the Parties do Make Measurements (PMM) as this affects the form of the probability to be calculated. 
                    % If k parties make measurements then there will be k products of variables.            
                    pmm = find(mvalues);
                    numpmm=length(pmm);
                    
                    cellsz = cellfun(@size, obj.coefflist(i2,:),'uni',false);
                    cellsize = cell2mat(cellsz);
                    lengthvar = cellsize(2);
                    
                    for i3 = 1:lengthvar
                        var = obj.coefflist{i2,1,1};
                        coeff = var(i3);
                        % If the coefficient is zero then just skip the calculation otherwise continue. 
                        if not(coeff) == 0           
                            % get the values of the outcomes observed from the values of i3 and the parties that make measurements.  
                            dvaluespmm = (dec2base(i3-1,obj.d,numpmm)+1)-'0';
                            dvalues = zeros(1,obj.n);
                            counter = 1;
                            %Add zeros back in to dvalues
                            for party = 1:obj.n
                                if ismember(party, pmm(:))
                                    dvalues(party) = dvaluespmm(counter);
                                    counter = counter + 1;
                                else
                                    dvalues(party) = 0;
                                end
                            end
                            % Now calculate the probability from the current location in the list i2,i3               
                            prob = calcprob(obj,pmm,mvalues,dvalues);                     
                            s = s + coeff*prob;
                        end
                    end
                end
                %tempsarray;
                % If smax is not yet defined then set the value of smax and the array of deterministic probabilities that give smax. 
                if obj.smax == 'x'
                  obj.smax = s;
                  obj.detprobsgivesmax(1,:) = obj.detprobvalues;
                  obj.detprobsrows = obj.detprobsrows + 1;
                else
                     % If s is greater than the current maximum then set the new value of smax and reset the array of deterministic 
                     % probabilities that give smax with this new value.
                     if (s > 1.0001*obj.smax)
                       obj.smax = s;
                       clear obj.detprobsgivesmax
                       obj.detprobsgivesmax = zeros(obj.d^(obj.n*obj.m),obj.numvars);
                       coder.varsize('obj.detprobsgivesmax',[obj.d^(obj.n*obj.m),obj.numvars]);
%                        obj.detprobsgivesmax = zeros(obj.maxdim,obj.numvars);
                       obj.detprobsgivesmax(1,:) = obj.detprobvalues;
                       obj.detprobsrows = 1;

                     % If s is the same as the current maximum, append the deterministic
                     % probabilities to the array of deterministic probabilities that give smax.
                     elseif ((0.9999*obj.smax <= s)  && (s <= 1.0001*obj.smax))
                       obj.detprobsgivesmax(obj.detprobsrows+1,:) = obj.detprobvalues;
                       obj.detprobsrows = obj.detprobsrows + 1;
                       
                     % Else If s is less than smax don't do anything, just continue to the next loop.          
                     else
                      
                     end    
                end  
          end
      end
      function index = getdetprobindex(obj,ni,di,mi)
          % GETDETPROBINDEX Calculates the index of the deterministic probability with D_ni(di,mi) within the array of deterministic probabilities.          
          index = (ni-1)*obj.d*obj.m+(di-1)*obj.m+mi; 
      end
      function prob = calcprob(obj,pmm,mvalues,dvalues)
      % CALCPROB Calculates the probability P(d1d2..dn|m1m2...mn)
      
           % The expression will take the form of k products of probabilities if there are k parties that make measurements. 
           %For each party that does make measurements calculate the corresponding product.
           % The NUMber of Parties that Make Measurements is the number of outcome variables di we have to loop over.              
           numpmm = length(pmm)
           prodterms = zeros(1,numpmm)
           for i1 = 1:numpmm
               party = pmm(i1)
               prodterms(i1) = obj.detprobvalues(getdetprobindex(obj,party,dvalues(party),mvalues(party)))
           end
           prob = prod(prodterms);
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
              prodterms = zeros(1,obj.n);
              for party = 1:obj.n
                  prodterms(party) = obj.detprobs(getdetprobindex(obj,party,obj.darray(party),obj.marray(party)));
              end
              curprodterm = prod(prodterms);
              obj.probdist(obj.probdistelements+1) = curprodterm;
              obj.probdistelements = obj.probdistelements + 1;
          end
      end
      function [dim] = calcdim(~,probdistsgivesmax)
      % CALCDIM Calculate the dimension of the bell inequality from the
      % probability distributions that give the classical bound.
      dim = rank(probdistsgivesmax); 
      end
   end
end