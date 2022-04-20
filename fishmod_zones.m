%function for fishing

function [catch_zones] = fishmod_zones ( N1, umat, spparams,dT)
    Fmortmat = umat{1,1}*spparams.q*dT; %fishing mortality matrix
	bio = N1{1,1} .* spparams.weightsmat; %biomass at year 1
	catch_zones = sum((bio .* Fmortmat .* (1 - exp(-Fmortmat - (spparams.M.*dT)))) ./ (Fmortmat + (spparams.M*dT))); % profits
end
