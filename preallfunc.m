%Function to preallocate and set abundance in first year

function[N1] = preallfunc(Ninit,spparams, umat, no_reefs, reefs_std)
	N1 = nan(spparams.nages,no_reefs); %initial matrix of fish age in all reefs
	agesrep = 0:(spparams.nages-1); %fish age range
	Fmortmat = umat{1,1}(:,:)*spparams.q;%fishing mortality rate, applied in all reefs
    mortrate = (exp(-agesrep .* (spparams.M+Fmortmat)'))';
       
    % N at year 1
	N1 = mortrate .* Ninit;  %calculate initial fmpa
    for i = 1:no_reefs
    N1(:,i) = N1(:,i)*reefs_std(i);
    end
    %initial abundance in reserve and unreserved in year 1 across all ages
end