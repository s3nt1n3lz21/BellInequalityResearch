classdef classicalBoundAndDimCalculator < handle
% CLASSICALBOUNDANDDIMCALCULATOR A class used to calculate the dimension and classical bound
% of a bell inequality.
   properties      
      % The number of parties in the scenario.
      noParties
      % The list of the total number of possible outcomes for each measurement of each party 
      maxNoMeasOutcomesList    
      % The list of probability coefficients that define the Bell operator.
      probCoeffList
      
      % A vector form of maxNoMeasOutcomesList
      maxNoMeasOutcomesVec       
      % The total number of local deterministic probabilities.
      totalNoLocalProbs
      % Holds numbers that determine which one of the local probabilities part of each constraint group is currently taking the value 1.
      indexArray
      % The maximum possible facet dimension of the Bell inequality (The spatial dimension).
      maxDim
      % An array holding the current values of the local deterministic probabilities.      
      localProbValues
      % The current maximum classical value of the Bell Inequality found so far.
      sMax
      % The array of local probabilities that give smax. Each row is one set of local probabilities that give sMax.      
      localProbsGiveSMax
      % How many rows full is the array of local probabilities (not the number of rows).     
      localProbsGiveSMaxRows 
      % An array used to hold the current probability distribution in calcProbDists
      probDist
      % The total number of constraints.
      noConstraints
      % A list of the maximum number of measurements that each party can make.
      maxMeasSettings
      % The length of each behaviour (probability vector)
      behaviourLength      
      % An array to hold a set of local probabilities as an array. Each row is the set of local probabilities for one of the parties. Used in CalcProbDists to calculate the behaviours. 
      localProbArray
      % A list used in calcProbDist to hold the current local probabilities of each party that will be multiplied together to calculate each element of the behaviour.
      localProbsToMultiply
      % A counter used in CalcProbDist to keep track of what element of the behaviour is being calculated.
      behaviourElementCounter
   end
   methods       
      function obj = classicalBoundAndDimCalculator(maxNoMeasOutcomesList,probCoeffList)
         % CBANDDIMCALC Constructor function to initialise the properties.
         obj.noParties = size(maxNoMeasOutcomesList,2);
         obj.maxNoMeasOutcomesList = maxNoMeasOutcomesList;
         obj.probCoeffList = probCoeffList;
         
         % Initialisation of temporary variables.
         % A variable used to calculate the spatial dimension of Bell inequalities of this scenario.
         prodVals = zeros(1,size(maxNoMeasOutcomesList,1));
         % A variable used to calculate the total number of local probabilities.
         totalNoLocalProbsTemp = 0;
         % An array used to calculate the behaviour length.
         noPossibleOutcomesArray = zeros(1,obj.noParties);       
         % Used to calculate the total number of possible measurements.
         noConstraintsTemp = 0;
         % Used to store the maximum number of measurements each party can make.
         maxMeasSettingsTemp = zeros(1,obj.noParties);
         
         % Loop over each party's list of outcomes to calculate and initialise some of the properties.
         for party = 1:size(maxNoMeasOutcomesList,2)
            currentMaxNoMeasOutcomesList = maxNoMeasOutcomesList{1,party,1};

            % Calculate the total number of local probabilities/outcomes
            totalNoLocalProbsTemp = totalNoLocalProbsTemp + sum(currentMaxNoMeasOutcomesList);
            % Used to calculate the spatial dimension.
            prodValue = 0;
            
            noPossibleOutcomesArray(party) = sum(currentMaxNoMeasOutcomesList);
            
            % Calculate the number of possible measurements for the current party (The length of the outcomes list)
            cellSz = cellfun(@size, maxNoMeasOutcomesList(:,party),'uni',false);
            cellSize = cell2mat(cellSz);
            currPartyTotalNoMeas = cellSize(1);
            
            %Store the maximum number of measurements that each party can make.
            maxMeasSettingsTemp(party) = currPartyTotalNoMeas;
                        
            % Calculate the total number of possible measurements.
            noConstraintsTemp = noConstraintsTemp + currPartyTotalNoMeas;
            
            % Create another list to hold the maximum number of measurement outcomes for each measurement of each party but as a 1-dimensional vector.
            if party == 1
                maxNoMeasOutcomesVecTemp = currentMaxNoMeasOutcomesList';
            else
                maxNoMeasOutcomesVecTemp = horzcat(maxNoMeasOutcomesVecTemp,currentMaxNoMeasOutcomesList');           
            end           
            
            % Calculate the spatial dimension
            for i2 = 1:currPartyTotalNoMeas
                prodValue = prodValue + currentMaxNoMeasOutcomesList(i2);
            end
            prodValue = prodValue - currPartyTotalNoMeas + 1;
            % Used to calculate the spatial dimension.
            prodVals(party) = prodValue;
                      
         end
         
         obj.maxDim = prod(prodVals) - 2;
         obj.behaviourLength = prod(noPossibleOutcomesArray);
         obj.maxNoMeasOutcomesVec = maxNoMeasOutcomesVecTemp;
         obj.totalNoLocalProbs = totalNoLocalProbsTemp;
         obj.localProbValues = zeros(1,obj.totalNoLocalProbs); 
         obj.noConstraints = noConstraintsTemp;
         obj.maxMeasSettings = maxMeasSettingsTemp;

         %Make detprobsgivesmax dynamic but specify an upper bound. The spatial dimension of the Bell inequality is an upper estimate of the maximum number of ways of getting smax, and so the memory allocated.
         obj.localProbsGiveSMax = zeros(obj.maxDim,obj.totalNoLocalProbs);
         coder.varsize('obj.localProbsGiveSMax',[obj.maxDim,obj.totalNoLocalProbs]);
         
         obj.localProbsGiveSMaxRows = 0;
         obj.sMax = 'x';
         obj.indexArray = zeros(1,obj.noConstraints);
         obj.localProbArray = {};
         
      end
      
      function [dim,sMax] = calc(obj)
      % CALC Start the calculation of the dimension and classical bound.
      
          % Calculate the maximum classical bound and local probabilities that give sMax.
          loopExtremalBehaviours(obj,obj.noConstraints);
          sMax = obj.sMax;
          
          % Calculate the behaviours from the local probabilities that give sMax.
          probDistsGiveSMax = calcProbDists(obj);
          
          % Calculate the facet dimension of the inequality.
          dim = calcDim(obj,probDistsGiveSMax);    
      end
      
      function loopExtremalBehaviours(obj,indexNosToLoop)
      % LOOPEXTREMALBEHAVIOURS The main algorithm to loop over all the possible extremal behaviours and calculate the Bell value.
      
          % The values of the index variables determine the values of the
          % local probabilities. If there is a still a index variable to loop over then loop
          % over its possible values.
          if indexNosToLoop >= 1
               maxIndexNoValue = obj.maxNoMeasOutcomesVec(obj.noConstraints-indexNosToLoop+1);
               for indexNo = 1:maxIndexNoValue
                  obj.indexArray(obj.noConstraints-indexNosToLoop+1) = indexNo;
                  % Then loop over the rest by calling the function again, with one less variable to loop over.    
                  loopExtremalBehaviours(obj,indexNosToLoop-1);
               end
               
          % Now calculate the Bell value for the current values of the local probabilities.
          else
                % Set the values of the local probabilities from the current values of the index variables.           
                obj.localProbValues = zeros(1,obj.totalNoLocalProbs);
                indexArrayIndex = 0;
                for party = 1:obj.noParties
                    maxMeasSetting = obj.maxMeasSettings(party);
                    for measSetting = 1:maxMeasSetting
                        indexArrayIndex = indexArrayIndex + 1;
                        indexNo = obj.indexArray(indexArrayIndex);
                        obj.localProbValues(getLocalProbIndex(obj,party,indexNo,measSetting)) = 1;
                    end
                end

                % Initialise the Bell value.              
                s = 0;
                             
                % For each term in the probability coefficient list calculate the corresponding contribution to the Bell value. e.g coeff*prob.
                              
                currMeasSettings = cat(2,zeros(1,obj.noParties-1),[-1]);
                for column = 1:size(obj.probCoeffList,2)
                    % Get the correct measurement settings (e.g m1 and m2) from how far into the coefficient list we are. 
                    % Calculate the measurement settings by keep adding one to measurement setting numbers. 

                    for party = obj.noParties:-1:1
                        % If the current measurement setting number is maximum set it to zero and carry it over to thenext parties measurement number.
                        if currMeasSettings(party) == obj.maxMeasSettings(party)
                            currMeasSettings(party) = 0;
                        else
                        % Add one to the measurement setting.    
                            currMeasSettings(party) = currMeasSettings(party) + 1;
                            break
                        end
                    end
                    
                    % Find which of the Parties do Make Measurements (PMM) as this affects the form of the probability to be calculated. 
                    % If k parties make measurements then there will be k products of local probabilities.       
                    partiesMakingMeas = find(currMeasSettings);
                    noPartiesMakingMeas=length(partiesMakingMeas);
                    
                    % Calculate the number of probability terms to be calculated.
                    cellSz = cellfun(@size, obj.probCoeffList(:,column),'uni',false);
                    cellSize = cell2mat(cellSz);
                    noProbTerms = cellSize(1);
                    
                    % Calculate the maximum number of outcomes for each party that is making a measurement. This is used to find which probabilities are going to be multiplied.               
                    outcomeNosOfPartiesMakingMeas = cat(2,ones(1,noPartiesMakingMeas-1),[0]);
                    maxOutcomeNosOfPartiesMakingMeas = zeros(1,noPartiesMakingMeas);
                    for index = 1:noPartiesMakingMeas
                        party = partiesMakingMeas(index);  
                        currPartysMaxNoMeasOutcomesList = obj.maxNoMeasOutcomesList{1,party,1};
                        maxOutcomeNo = currPartysMaxNoMeasOutcomesList(currMeasSettings(party));
                        maxOutcomeNosOfPartiesMakingMeas(index) = maxOutcomeNo;
                    end
                    
                    for row = 1:noProbTerms
                        % Calculate which of the local probabilities are going to be multiplied from how far into this list of probabilities we are, for these measurement settings.
                        
                        % Calculate the coeffient of the probability term.
                        var = obj.probCoeffList{1,column,1};
                        coeff = var(row);
                        
                        % If the coefficient is zero then just skip the calculation otherwise continue. 
                        if not(coeff) == 0
                            
                            % Calculate which of the local probabilities are going to be multiplied by keep adding one to the list of outcome numbers.                            
                            for index = noPartiesMakingMeas:-1:1
                                if outcomeNosOfPartiesMakingMeas(index) == maxOutcomeNosOfPartiesMakingMeas(index)
                                    outcomeNosOfPartiesMakingMeas(index) = 1;
                                else
                                    outcomeNosOfPartiesMakingMeas(index) = outcomeNosOfPartiesMakingMeas(index) + 1;
                                    break
                                end
                            end
                            
                            % Calculate the final list of outcome numbers by adding zeros back in for the parties that don't make measurements.
                            outcomeNos = zeros(1,obj.noParties);
                            partiesMakingMeasCounter = 1;
                            for party = 1:obj.noParties
                                if ismember(party, partiesMakingMeas(:))
                                    outcomeNos(party) = outcomeNosOfPartiesMakingMeas(partiesMakingMeasCounter);
                                    partiesMakingMeasCounter = partiesMakingMeasCounter + 1;
                                else
                                    outcomeNos(party) = 0;
                                end
                            end
                            
                            % Now calculate the probability from the current measurement settings and outcome numbers.           
                            prob = calcProb(obj,partiesMakingMeas,currMeasSettings,outcomeNos);
                            
                            % Calculate the contribution to the Bell value.
                            s = s + coeff*prob;
                        end
                    end
                end
                % If the maximum Bell value is not yet defined then set its value to the result of the first calculation. Store the set of local probabilities that give the current maximum Bell value. 
                if obj.sMax == 'x'
                  obj.sMax = s;
                  obj.localProbsGiveSMax(1,:) = obj.localProbValues;
                  obj.localProbsGiveSMaxRows = obj.localProbsGiveSMaxRows + 1;
                else
                    
                     % If the current Bell value is greater than the current maximum then set the new value of the maximum and clear the list of local probabilities that give the maximum.
                     if (s > 1.0001*obj.sMax)
                       obj.sMax = s;
                       clear obj.localProbsGiveSMax
                       obj.localProbsGiveSMax = zeros(obj.noConstraints,obj.totalNoLocalProbs);
                       coder.varsize('obj.localProbsGiveSMax',[obj.noConstraints,obj.totalNoLocalProbs]);
                       
                       obj.localProbsGiveSMax(1,:) = obj.localProbValues;
                       obj.localProbsGiveSMaxRows = 1;

                     % If the current Bell value is the same as the current maximum add the local probability to the list of probabilities that give the maximum.
                     elseif ((0.9999*obj.sMax <= s)  && (s <= 1.0001*obj.sMax))
                       obj.localProbsGiveSMax(obj.localProbsGiveSMaxRows+1,:) = obj.localProbValues;
                       obj.localProbsGiveSMaxRows = obj.localProbsGiveSMaxRows + 1;

                     % If the current Bell value is less than the maximum then just continue to the next extremal point.         
                     else
                      
                     end    
                end  
          end
      end
      
      function index = getLocalProbIndex(obj,partyNo,outcomeNo,measSettingNo)
      % GETLOCALPROBINDEX Calculates the index of the local probability of party (partyNo) doing measurement (measSetting) and getting outcome (outcomeNo) e.g D_ni(di,mi) within the list of local probabilities
          
          % The index is found by summing the outcome numbers.
          index = 0;
          for party = 1:partyNo-1
            currPartysMaxNoMeasOutcomesList = obj.maxNoMeasOutcomesList{1,party,1};
            index = index + sum(currPartysMaxNoMeasOutcomesList);
          end  
          
          currPartysMaxNoMeasOutcomesList = obj.maxNoMeasOutcomesList{1,partyNo,1};
            
          for measSetting = 1:measSettingNo-1
            index = index + currPartysMaxNoMeasOutcomesList(measSetting);
          end
          
          index = index + outcomeNo;
      end
      
      function prob = calcProb(obj,partiesMakingMeas,currMeasSettings,outcomeNos)
      % CALCPROB Calculates the probability that the parties get the given outcomes for the measurements they make P(d1d2..dk|m1m2...mk)
      
           % The expression will take the form of k products of probabilities if there are k parties that make measurements. 
           % For each party that does make a measurement calculate the corresponding local probability they get their particular outcome.           
           noPartiesMakingMeas = length(partiesMakingMeas);
           localProbsBeingMultiplied = zeros(1,noPartiesMakingMeas);
           for counter = 1:noPartiesMakingMeas
               party = partiesMakingMeas(counter);
               localProbsBeingMultiplied(counter) = obj.localProbValues(getLocalProbIndex(obj,party,outcomeNos(party),currMeasSettings(party)));
           end
           
           % Multiply the local probabilities together.
           prob = prod(localProbsBeingMultiplied);
      end
      
      function [probDistsGiveSMax] = calcProbDists(obj)
      % CALCPROBDISTS Calculate the behaviours P(d1d2..dn|m1m2...mn) that give the classical bound from the array of local probabilities that give the classical bound.
          
          % For each set of local deterministic probabilities that give the bound calculate the behaviour.
          probDistsGiveSMax = zeros(obj.behaviourLength,obj.localProbsGiveSMaxRows);
          for row = 1:obj.localProbsGiveSMaxRows
              obj.probDist = zeros(1,obj.behaviourLength);
              localDetProbs = obj.localProbsGiveSMax(row,:);
         
              % Calculate the behaviour by looping over all the possible products of local probabilities P=D1D2...Dn for each D1,D2,...,Dn
              % Do this by splitting up the list of local probabilities into lists of local probabilities for each party.
              sliceOfLocalProbs = localDetProbs;
              for party = 1:obj.noParties   
                currPartysMaxNoMeasOutcomesList = obj.maxNoMeasOutcomesList{1,party,1};
                currPartysLocalProbListSize = sum(currPartysMaxNoMeasOutcomesList);
                obj.localProbArray(party,:) = {sliceOfLocalProbs(1:currPartysLocalProbListSize)};
                sliceOfLocalProbs = sliceOfLocalProbs(currPartysLocalProbListSize+1:length(sliceOfLocalProbs));
              end
              
              % Calculate the current behaviour/probability distribution.
              obj.localProbsToMultiply = zeros(1,obj.noParties);
              obj.behaviourElementCounter = 1;
              calcProbDist(obj,obj.noParties);
              % Store the probability distributions for calculation of the
              % dimension. Note we want the probability distributions to be
              % the columns of the matrix, not the rows.
              column = row;
              probDistsGiveSMax(:,column) = obj.probDist;
          end
      end
      
      function calcProbDist(obj,partiesToLoop)
      % CALCPROBDIST Calculate the probability distribution P(d1d2..dn|m1m2...mn) from the local deterministic probabilities.
      
          % Loop over all the possible products of the local probabilities of each party through recursion.
          if partiesToLoop >= 1
              for loopcounter = 1:length(obj.localProbArray{obj.noParties-partiesToLoop+1,1,1})
                  currPartysLocalProbList = obj.localProbArray{obj.noParties-partiesToLoop+1,1,1};
                  obj.localProbsToMultiply(obj.noParties-partiesToLoop+1) = currPartysLocalProbList(loopcounter);
                  calcProbDist(obj,partiesToLoop-1);       
              end
          else
              % Calculate the product of the current local probabilities being multiplied.        
              obj.probDist(obj.behaviourElementCounter) = prod(obj.localProbsToMultiply);
              obj.behaviourElementCounter = obj.behaviourElementCounter + 1;
          end
      end
      
      function [dim] = calcDim(~,probDistsGiveSMax)
      % CALCDIM Calculate the dimension of the bell inequality from the probability distributions that give the classical bound.
          dim = rank(probDistsGiveSMax)-1;
      end
      
   end
end

%TODO LIST:
%Do validation on size and shape of probCoeffList based on maxNoMeasOutcomesList
%Update ReadMe