classdef cbanddimcalc < handle
% CBANDDIMCALC A class used to calculate the dimension and classical bound
% of a bell inequality.
   properties      
      % The number of parties n, measurement outcomes d and measurement settings m. Here they are assumed to all be the same.      
      n
      %d
      %m
      % The number of single party deterministic probabilities.
      numVars
      % The array to hold the di values of each constraint equation that are 1.
      indexArray
      % The maximum possible dimension of the bell inequality is an estimate of
      % the maximum number of ways of getting smax, and so the memory allocated.
      maxDim
      % The correlator coefficient list.  
      coeffList
      % The array of the current values of the single party deterministic probabilities.      
      detProbValues
      % The current maximum found so far.
      sMax
      % The array of deterministic probabilities that give smax.      
      detProbsGiveSMax
      % How many rows full is the array of deterministic probabilities (not the number of rows).     
      detProbsRows
      % An array used in calcprobdists to hold the current values of the possible
      % outcomes.
      dArray
      % An array used in calcprobdists to hold the current values of the
      % measurent settings.
      mArray
      % An array used to hold the current probability distribution in
      % calcProbDists
      probDist
      % How full is the vector probdist
      probDistElements
      % The list of the total number of possible outcomes for each measurement of each party 
      dList
      % A vector form of the dlist
      dVecList
      % The total number of possible measurements
      totalNumMeas
      % The total length of each probability vector
      maximumDetProbs
      % The local deterministic probabilities as an array, each row is the
      % local deterministic probabilities of one of the parties
      dProbArray
      % A list used in calcProbDist to hold the d values that are looped over
      currDVals
      % A counter used in CalcProbDist for the index of the current probability term in
      % the probability distribution
      probCounter
   end
   methods       
      function obj = cbanddimcalc(n,dList,coeffList)
         % CBANDDIMCALC Constructor function to the initialise properties.
         obj.n = n;
         obj.dList = dList;
         % Calculate the spatial dimension of bell inequalities of this scenario
         prodVals = zeros(1,size(dList,1));
         % The total number of deterministic probabilities
         numVars = 0;
         % An array used to calculate the maximum possible deterministic probabilities "maximumdetprobs"
         tempVarArray = zeros(1,n);       
         % The total number of possible measurements to be used later for "indexarray"
         totalNumMeas = 0;
         for i1 = 1:size(dList,1)
            var = dList{i1,1,1};

            % Calculate the total number of deterministic probabilities
            numVars = numVars + sum(var);
            
            prodValue = 0;
            tempVarArray(i1) = sum(var);
            
            cellSz = cellfun(@size, dList(i1,:),'uni',false);
            cellSize = cell2mat(cellSz);
            lengthVar = cellSize(2);
            
            % Calculate the total number of possible measurements to be used later for "indexarray"
            totalNumMeas = totalNumMeas + lengthVar;
            
            if i1 == 1
                dVecList = var;
            else
                dVecList = horzcat(dVecList,var);           
            end           
            
            for i2 = 1:lengthVar
                prodValue = prodValue + var(i2);
            end
            prodValue = prodValue - lengthVar + 1;
            prodVals(i1) = prodValue;
                      
         end                             
         obj.maxDim = prod(prodVals) - 1;
         obj.maximumDetProbs = prod(tempVarArray);
         obj.dVecList = dVecList;
         %obj.d = d;
         %obj.m = m;
         
         obj.numVars = numVars;
         obj.totalNumMeas = totalNumMeas;
         obj.coeffList = coeffList;
         obj.detProbValues = zeros(1,obj.numVars);
         %Make detprobsgivesmax dynamic but specify an upper bound.
         obj.detProbsGiveSMax = zeros(obj.maxDim,obj.numVars);
         coder.varsize('obj.detProbsGiveSmax',[obj.maxDim,obj.numVars]);
         obj.detProbsRows = 0;
         obj.sMax = 'x';
         obj.indexArray = zeros(1,obj.maxDim);
         obj.dProbArray = {};
         %listsize = size(coefflist,1);
         %if not(listsize == (m+1)^n)
         %   fprintf("Error, the dimension of the correlator coefficient list does not match the scenario\n The dimension is %.0f when it should be %.0f",size(corrcoefflist,2),(m+1)^n)
         %end
      end
      function [dim,sMax] = calc(obj)
      % CALC Calculate the dimension and classical bound.
       loopDetProbs(obj,obj.totalNumMeas);
       % Calculate the probability distribution vectors from determinsitic
       % probability array.
       sMax = obj.sMax;
       probDistsGiveSMax = calcProbDists(obj);
       % Calculate the dimension.
       dim = calcDim(obj,probDistsGiveSMax);    
      end
      function loopDetProbs(obj,varsToLoop)
      % LOOPDETPROBS The main algorithm to loop over the possible values of the
      % deterministic probabilities (0 or 1).
      
          % If there is a still a deterministic probability to loop over then loop
          % over its possible values.
          if varsToLoop >= 1
               maxDValue = obj.dVecList(obj.totalNumMeas-varsToLoop+1); %%%%%%%%%%%%Does the order matter?
               for value = 1:maxDValue
                  obj.indexArray(obj.totalNumMeas-varsToLoop+1) = value;
                  % Then loop over the rest by calling the function again, with one less variable to loop over.    
                  loopDetProbs(obj,varsToLoop-1);
               end
               
          % If there are no more deterministic probabilities to loop over then
          % calculate the Bell Inequality correlator for the current values of the
          % deterministic probabilities.
          else
                % Initialise the value of the expression to be 0 and get the values of the variables from the indexes.         
                obj.detProbValues = zeros(1,obj.numVars);
                dIndex = 0;
                for party = 1:obj.n
                    cellSz = cellfun(@size, obj.dList(party,:),'uni',false);
                    cellSize = cell2mat(cellSz);
                    lengthVar = cellSize(2);
                    for mi = 1:lengthVar
                        dIndex = dIndex + 1;
                        di = obj.indexArray(dIndex);
                        obj.detProbValues(getDetProbIndex(obj,party,di,mi)) = 1;
                    end
                end

                s = 0;
                             
                % For each term in the coefficient list calculate the corresponding contribution from that probability. e.g coeff*prob
                % Get the correct measurement settings (e.g m1 and m2) from how far into the coefficient list we are by firstcalculating the maximum m values.
                mValues = cat(2,zeros(1,obj.n-1),[-1]);
                maxMValues = zeros(1,obj.n);
                for dIndex = 1:obj.n
                    cellSz = cellfun(@size, obj.dList(dIndex,:),'uni',false);
                    cellSize = cell2mat(cellSz);
                    maxMValue = cellSize(2);
                    maxMValues(dIndex) = maxMValue;
                end    
                for i2 = 1:size(obj.coeffList,1)
                    % Get the correct measurement settings (e.g m1 and m2) from how far into the coefficient list we are
                    %mvalues = dec2base(i2-1,obj.m+1,obj.n)-'0';

                    for index = obj.n:-1:1
                        if mValues(index) == maxMValues(index)
                            mValues(index) = 0;
                        else
                            mValues(index) = mValues(index) + 1;
                            break
                        end
                    end
                    
                    % Find which of the Parties do Make Measurements (PMM) as this affects the form of the probability to be calculated. 
                    % If k parties make measurements then there will be k products of variables.            
                    partiesMakeMeas = find(mValues);
                    numPartiesMakeMeas=length(partiesMakeMeas);
                    
                    cellSz = cellfun(@size, obj.coeffList(i2,:),'uni',false);
                    cellSize = cell2mat(cellSz);
                    lengthVar = cellSize(2);
                    
                    dValuesPMM = cat(2,ones(1,numPartiesMakeMeas-1),[0]);
                    maxDValuesPMM = zeros(1,numPartiesMakeMeas);
                    for dIndex = 1:numPartiesMakeMeas
                        party = partiesMakeMeas(dIndex);  
                        dvec = obj.dList{party,1,1};
                        maxDValue = dvec(mValues(party));
                        maxDValuesPMM(dIndex) = maxDValue;
                    end
                    for i3 = 1:lengthVar
                        var = obj.coeffList{i2,1,1};
                        coeff = var(i3);
                        % If the coefficient is zero then just skip the calculation otherwise continue. 
                        if not(coeff) == 0
                            for index = numPartiesMakeMeas:-1:1
                                if dValuesPMM(index) == maxDValuesPMM(index)
                                    dValuesPMM(index) = 1;
                                else
                                    dValuesPMM(index) = dValuesPMM(index) + 1;
                                    break
                                end
                            end
                            % get the values of the outcomes observed from the values of i3 and the parties that make measurements.  
                            dValues = zeros(1,obj.n);
                            counter = 1;
                            %Add zeros back in to dvalues
                            for party = 1:obj.n
                                if ismember(party, partiesMakeMeas(:))
                                    dValues(party) = dValuesPMM(counter);
                                    counter = counter + 1;
                                else
                                    dValues(party) = 0;
                                end
                            end
                            %dvalues
                            % Now calculate the probability from the current location in the list i2,i3               
                            prob = calcProb(obj,partiesMakeMeas,mValues,dValues);                     
                            s = s + coeff*prob;
                        end
                    end
                end
                s;
                %tempsarray;
                % If smax is not yet defined then set the value of smax and the array of deterministic probabilities that give smax. 
                if obj.sMax == 'x'
                  obj.sMax = s;
                  obj.detProbsGiveSMax(1,:) = obj.detProbValues;
                  obj.detProbsRows = obj.detProbsRows + 1;
                else
                     % If s is greater than the current maximum then set the new value of smax and reset the array of deterministic 
                     % probabilities that give smax with this new value.
                     if (s > 1.0001*obj.sMax)
                       obj.sMax = s;
                       clear obj.detProbsGiveSMax
                       obj.detProbsGiveSMax = zeros(obj.totalNumMeas,obj.numVars);
                       coder.varsize('obj.detProbsGiveSmax',[obj.totalNumMeas,obj.numVars]);
                       
                       obj.detProbsGiveSMax(1,:) = obj.detProbValues;
                       obj.detProbsRows = 1;

                     % If s is the same as the current maximum, append the deterministic
                     % probabilities to the array of deterministic probabilities that give smax.
                     elseif ((0.9999*obj.sMax <= s)  && (s <= 1.0001*obj.sMax))
                       obj.detProbsGiveSMax(obj.detProbsRows+1,:) = obj.detProbValues;
                       obj.detProbsRows = obj.detProbsRows + 1;

                     % Else If s is less than smax don't do anything, just continue to the next loop.          
                     else
                      
                     end    
                end  
          end
      end
      function index = getDetProbIndex(obj,ni,di,mi)
          % GETDETPROBINDEX Calculates the index of the deterministic probability with D_ni(di,mi) within the array of deterministic probabilities.    
          index = 0;
          for i1 = 1:ni-1
            var = obj.dList{i1,1,1};
            index = index + sum(var);
          end  
            
          %cellsz = cellfun(@size, dlist(i1,:),'uni',false);
          %cellsize = cell2mat(cellsz);
          %lengthvar = cellsize(2);
          
          var = obj.dList{ni,1,1};
            
          for i2 = 1:mi-1
            index = index + var(i2);
          end
          
          index = index + di;
          
          %index = (ni-1)*obj.d*obj.m+(di-1)*obj.m+mi; %%%AS VECTOR ORDERS THEM BY N then D then M this is not correct?
      end
      function prob = calcProb(obj,pMM,mValues,dValues)
      % CALCPROB Calculates the probability P(d1d2..dn|m1m2...mn)
      
           % The expression will take the form of k products of probabilities if there are k parties that make measurements. 
           %For each party that does make measurements calculate the corresponding product.
           % The NUMber of Parties that Make Measurements is the number of outcome variables di we have to loop over.              
           numPMM = length(pMM);
           prodTerms = zeros(1,numPMM);
           for i1 = 1:numPMM
               party = pMM(i1);
               prodTerms(i1) = obj.detProbValues(getDetProbIndex(obj,party,dValues(party),mValues(party)));
           end
           prob = prod(prodTerms);
      end
      function [probDistsGiveSMax] = calcProbDists(obj)
      % CALCPROBDISTS Calculate the probability distribution vectors that give
      % the classical bound from the array of deterministic probabilities that give the classical bound.
          obj.mArray = zeros(1,obj.n);
          obj.dArray = zeros(1,obj.n);
          probDistsGiveSMax = zeros(obj.detProbsRows,obj.maximumDetProbs);
          %For each list of local deterministic probs that gives smax
          %calculate the behaviour
          
          %There is only going to be one index in the probability
          %distribution that is 1.
          for row = 1:obj.detProbsRows
              obj.probDist = zeros(1,obj.maximumDetProbs);
              localDetProbs = obj.detProbsGiveSMax(row,:);
              %Calculate all the possible products of the local probabilities
              %Split the dVecList up into lists of d values for each party
              sliceOfDetProbs = localDetProbs;
              for i1 = 1:obj.n   
                var = obj.dList{i1,1,1};
                dListSize = sum(var);
                obj.dProbArray(i1,:) = {sliceOfDetProbs(1:dListSize)};
                sliceOfDetProbs = sliceOfDetProbs(dListSize+1:length(sliceOfDetProbs));
              end
              
              obj.currDVals = zeros(1,obj.n);
              obj.probCounter = 1;
              calcProbDist(obj,obj.n);    
              probDistsGiveSMax(row,:) = obj.probDist;
          end
      end
      function calcProbDist(obj,partiesToLoop)
      % CALCPROBDIST Calculate the probability distribution from the local
      % deterministic probabilities.
      
      if partiesToLoop >= 1
          for loopcounter = 1:length(obj.dProbArray{obj.n-partiesToLoop+1,1,1})
              currPartyDList = obj.dProbArray{obj.n-partiesToLoop+1,1,1};
              obj.currDVals(obj.n-partiesToLoop+1) = currPartyDList(loopcounter);
              calcProbDist(obj,partiesToLoop-1);       
          end
      else
          obj.probDist(obj.probCounter) = prod(obj.currDVals);
          obj.probCounter = obj.probCounter + 1;
          %%%%%%%%%%%%%%%%Set each element of prob dist
      end
      end
      function [dim] = calcDim(~,probDistsGiveSMax)
      % CALCDIM Calculate the dimension of the bell inequality from the
      % probability distributions that give the classical bound.
      dim = rank(probDistsGiveSMax); 
      end
   end
end