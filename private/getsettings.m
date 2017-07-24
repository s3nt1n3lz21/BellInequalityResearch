function settings = getsettings(index,m,n)
% function to calculate what the measurements settings are for the current correlator
    % create an array to hold these values    
    marray = zeros(1,n);
    % in general the index of a particular group of settings is given by "index = 1 + m1*(m+1)^(n-1) + m2*(m+1)^(n-2) + m3*(m+1)^(n-3)+..." 
    remainder = (index-1)/((m+1)^(n-1));
    for i1 = 1:n
        while remainder >= 1-eps(1)
            marray(i1) = marray(i1) + 1;
            remainder = remainder - 1;
        end
        remainder = remainder*(m+1);
    end
    settings = marray;    
    