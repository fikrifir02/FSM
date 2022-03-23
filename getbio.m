%function to calculate biomass and  catch
function[Biom, Catch]= getbio (N1, spparams)
	Biom = sum(N1{1,1} .* spparams.weightsmat); %To calculate biomass
end 