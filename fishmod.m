%function for fishing

function [Ben] = fishmod ( N1, umat, spparams, reefs_std,dT)
	Fmortmat = umat{1,1}*spparams.q*dT; %fishing mortality matrix
	bio = N1{1,1} .* spparams.weightsmat; %biomass at year 1
	profs = sum((bio .* Fmortmat .* (1 - exp(-Fmortmat - (spparams.M.*dT)))) ./ (Fmortmat + (spparams.M*dT)) * spparams.RelVal); % profits
	costs = profs .* (spparams.omega./(sum(bio) ./ reefs_std')); %cost is used to calculate profits/Benefits
	Ben = sum(profs(isfinite(costs))-costs(isfinite(costs))); %Ben, benefits is used to calculate harvest
end
