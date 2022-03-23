%Function to preallocate and set abundance in first year

function[N1] = preallfunc(Ninit,spparams, A, umat)
	N1 = nan(spparams.nages,1); %initial matrix of fish age in nmpa, reserve and fmpa
	agesrep = 0:(spparams.nages-1); %fish age range
	Fmortmat = umat{1,1}(:,2)*spparams.q;%fishing mortality rate, applied in non reserve area
    mortrate = exp(-agesrep .* (spparams.M+Fmortmat)');
            if 0 <= A(1) + A(2) + A(3)<= 1 
	N1(:,3)= Ninit *  A(3)  * mortrate;  %calculate initial fmpa
	N1(:,2)= Ninit *  A(2) * mortrate; %calculate initial reserve
    N1(:,1)= Ninit * (1 - A(2)-A(3)) * mortrate;
    N1; %initial abundance in reserve and unreserved in year 1 across all ages
            else
            end
end