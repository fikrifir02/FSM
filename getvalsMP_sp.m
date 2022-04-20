%function to collect long term catch (tmax) for each runs

function[catch_opt] = getvalsMP_sp(nvals,no_reefs, reefs_status,reefs_conn, reefs_std, Etot, index_sp,spparams, Ninit, tmax, startmpa,Frem)
	catch_sp = nan(1,nvals);
	for irun = 1:nvals
		umat = arrayfun(@(index_sp) createumat(index_sp, spparams.afishind, Etot),index_sp, 'UniformOutput',false);
		xres = timeloopfunc(no_reefs, reefs_status,reefs_conn, reefs_std, umat, index_sp, Ninit,spparams, tmax, startmpa,Frem);
		catch_sp = sum(xres.Catches(:,:),2);
        catch_opt(irun) = catch_sp(tmax);
        catch_opt;
    end
end