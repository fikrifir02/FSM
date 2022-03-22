%select species parameters from spparams object
function [sp_spparams] = select_spparams (spparams, ispp)
spparams.nages = spparams.nages(ispp);
spparams.q = spparams.q(ispp);
spparams.M = spparams.M (ispp);
spparams.matind = spparams.matind{1,ispp};
spparams.weights = spparams.weights{1,ispp};
spparams.weightsmat = spparams.weightsmat{1,ispp};
spparams.fecparam = spparams.fecparam (ispp);
spparams.CR = spparams.CR(ispp);
spparams.Rmax = spparams.Rmax(ispp);
spparams.afishind = spparams.afishind{1,ispp};
spparams.RelVal = spparams.RelVal(ispp);
spparams.alpha = spparams.alpha(ispp);
spparams.beta = spparams.beta(ispp);
spparams.Lzero = spparams.Lzero(ispp);
spparams.mrate = spparams.mrate(ispp);
spparams.omega = spparams.omega (ispp);
spparams.minbio_prof = spparams.minbio_prof(ispp);
spparams.nmonths = spparams.nmonths;
spparams.dT = spparams.dT;
sp_spparams = spparams;
end