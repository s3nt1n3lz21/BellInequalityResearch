classdef classicalboundcalculator
   properties
      detprobvalues
      detprobsarray
   end
   methods
      function r = calcclassicalbound(obj)
         r = round([obj.Value],2);
      end
      function r = loopdetprobs(obj,n)
         r = [obj.Value] * n;
      end
   end
end