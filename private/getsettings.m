function settings = getsettings(index,m,n)
% function to calculate what the measurements settings are for the current correlator
    % create an array to hold these values    
    marray = zeros(1,n);
    % in general the index of a particular group of settings is given by "index = 1 + m1*(m+1)^(n-1) + m2*(m+1)^(n-2) + m3*(m+1)^(n-3)+..."    
    marray(1) = idivide(int32(index-1),(m+1)^(n-1));
    remainder = (index-1)/((m+1)^(n-1))-marray(1);
    for i1 = 2:n
        marray(i1) = round(remainder*(m+1));
        remainder = remainder*(m+1) - marray(i1);
    end
    settings = marray;       