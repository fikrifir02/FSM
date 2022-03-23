%function to collect long term catch (tmax) for each runs

function[catch_opt] = getvalsMP_sp(nvals,Aval,Etot, index_sp,spparams,Ainit, Ninit, tmax, startmpa)
	catch_sp = nan(1,nvals);
	for irun = 1:nvals
		umat = arrayfun(@(index_sp) createumat(index_sp, spparams.afishind, Etot),index_sp, 'UniformOutput',false);
		xres = timeloopfunc(Aval, umat, index_sp, Ninit,spparams, Ainit, tmax, startmpa);
		[catch_species catch_sp] = catch_per_species(xres);
        catch_opt(irun) = catch_sp(tmax);
        catch_opt;
    end
end