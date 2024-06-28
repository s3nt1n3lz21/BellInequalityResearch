classdef extremalPointLooper < handle
    properties 
        % The total number of parties part of the scenario
        noParties
        % The number of constraintn equations
        noConstraints
        % Holds numbers that determine which one of the local probabilities part of each constraint group is currently taking the value 1.
        indexArray
        % Holds the current values of the local probabilities
        localProbValues
        % Holds the maximum number of measurement settings for each party
        maxMeasSettings
        % An array that holds the number of maximum number of measurement
        % outcomes
        maxNoMeasOutcomesList
        % The number of elements of the Behaviours
        behaviourLength
        % An array to hold a set of local probabilities as an array. Each row is the set of local probabilities for one of the parties. Used in CalcProbDist to calculate 
        localProbArray
        % A list used in calcProbDist to hold the current local probabilities of each party that will be multiplied together to calculate each element of the behaviour.
        localProbsToMultiply
        % A counter used in CalcProbDist to keep track of what element of the behaviour is being calculated.
        behaviourElementCounter
        % An array used to hold the current probability distribution in calcProbDist
        probDist
        % The matrix A which has columns that are the behaviours of the extremal points
        matrixA
        % A counter that is used to signify how full the matrixA is
        matrixAColumnCounter
        % A vector form of maxNoMeasOutcomesList
        maxNoMeasOutcomesVec  
        % The total number of local deterministic probabilities.
        totalNoLocalProbs
    end
    methods
        function obj = extremalPointLooper(maxNoMeasOutcomesList)
            obj.maxNoMeasOutcomesList = maxNoMeasOutcomesList;
            obj.noParties = size(maxNoMeasOutcomesList,2);
            
            % Calculate the number of possible measurements for each party (The length of the outcomes list)
            cellSz = cellfun(@size, maxNoMeasOutcomesList(:,1),'uni',false);
            cellSize = cell2mat(cellSz);
            m = cellSize(1);
            
            % Calculate d
            dList = maxNoMeasOutcomesList{1,1,1};
            d = dList(1);
            
            obj.noConstraints = obj.noParties*m;
            obj.indexArray = zeros(1,obj.noConstraints);
            obj.maxMeasSettings = ones(1,obj.noParties)*m;
            obj.localProbArray = {};
            totalNoLocalProbsTemp = 0;
            % An array used to calculate the behaviour length.
            noPossibleOutcomesArray = zeros(1,obj.noParties); 
            
            % Loop over each party's list of outcomes to calculate and initialise some of the properties.
            for party = 1:size(maxNoMeasOutcomesList,2)
               currentMaxNoMeasOutcomesList = maxNoMeasOutcomesList{1,party,1};
               
               % Calculate the total number of local probabilities/outcomes
               totalNoLocalProbsTemp = totalNoLocalProbsTemp + sum(currentMaxNoMeasOutcomesList);
               
               % Calculate the length of the behaviours
               noPossibleOutcomesArray(party) = sum(currentMaxNoMeasOutcomesList);
               
               % Create another list to hold the maximum number of measurement outcomes for each measurement of each party but as a 1-dimensional vector.
               if party == 1
                   maxNoMeasOutcomesVecTemp = currentMaxNoMeasOutcomesList';
               else
                   maxNoMeasOutcomesVecTemp = horzcat(maxNoMeasOutcomesVecTemp,currentMaxNoMeasOutcomesList');           
               end 
            end
            
            obj.behaviourLength = prod(noPossibleOutcomesArray);
            obj.matrixA = zeros(obj.behaviourLength,d^(obj.noParties*m));
            obj.matrixAColumnCounter = 1;
            obj.totalNoLocalProbs = totalNoLocalProbsTemp;
            obj.maxNoMeasOutcomesVec = maxNoMeasOutcomesVecTemp;
        end
        function [A] = calc(obj)
        % CALC Start the calculation of the extremal points and matrix A.

          % Start looping over the extremal Behaviours and calculate A
          loopExtremalBehaviours(obj,obj.noConstraints);

          % Return the matrix.
          A = obj.matrixA;   
        end
        function loopExtremalBehaviours(obj,indexNosToLoop)
        % LOOPEXTREMALBEHAVIOURS A function to loop over all the possible extremal behaviours and calculate the Bell value

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
          % If no more index variables to loop over then calculate the
          % behaviour of this extremal point.
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
                
                % Calculate the behaviour (probability distribution vector) from the values of the local probabilities
                calcProbDist(obj);
                %obj.probDist
                obj.matrixA(:,obj.matrixAColumnCounter) = obj.probDist;
                obj.matrixAColumnCounter = obj.matrixAColumnCounter + 1;
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
        function calcProbDist(obj)
        % CALCPROBDIST Calculate the behaviour P(d1d2..dn|m1m2...mn) of the
        % extremal point from the local probabilities.
          
          % Initialise the probability distribution vector
          obj.probDist = zeros(1,obj.behaviourLength)';

          % Calculate the behaviour by looping over all the possible products of local probabilities P=D1D2...Dn for each D1,D2,...,Dn
          % Do this by splitting up the list of local probabilities into lists of local probabilities for each party.
          sliceOfLocalProbs = obj.localProbValues;
          for party = 1:obj.noParties   
            currPartysMaxNoMeasOutcomesList = obj.maxNoMeasOutcomesList{1,party,1};
            currPartysLocalProbListSize = sum(currPartysMaxNoMeasOutcomesList);
            obj.localProbArray(party,:) = {sliceOfLocalProbs(1:currPartysLocalProbListSize)};
            sliceOfLocalProbs = sliceOfLocalProbs(currPartysLocalProbListSize+1:length(sliceOfLocalProbs));
          end

          % Calculate the current behaviour/probability distribution.
          obj.localProbsToMultiply = zeros(1,obj.noParties);
          obj.behaviourElementCounter = 1;
          calcProbDistElements(obj,obj.noParties);

        end
        function calcProbDistElements(obj,partiesToLoop)
        % CALCPROBDISTELEMENTS Calculate the probability distribution P(d1d2..dn|m1m2...mn) from the local deterministic probabilities.
      
          % Loop over all the possible products of the local probabilities of each party through recursion.
          if partiesToLoop >= 1
              for loopcounter = 1:length(obj.localProbArray{obj.noParties-partiesToLoop+1,1,1})
                  currPartysLocalProbList = obj.localProbArray{obj.noParties-partiesToLoop+1,1,1};
                  obj.localProbsToMultiply(obj.noParties-partiesToLoop+1) = currPartysLocalProbList(loopcounter);
                  calcProbDistElements(obj,partiesToLoop-1);       
              end
          else
              % Calculate the product of the current local probabilities being multiplied.        
              obj.probDist(obj.behaviourElementCounter) = prod(obj.localProbsToMultiply);
              obj.behaviourElementCounter = obj.behaviourElementCounter + 1;
          end
      end
    end
end