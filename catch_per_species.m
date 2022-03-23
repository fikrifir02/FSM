%function to accumulate catches for each species in all zones

function [catch_species catch_sp] = catch_per_species(xres)
            catch_species = struct();
            catch_species.species = [];
            catch_species.sum_species = [];
            catch_sp = [];

            catch_species.species = xres.Catches(:,1+(3*(1-1)):3*1);
            catch_species.sum_species = sum(catch_species.species,2);
            catch_sp = [catch_sp catch_species.sum_species];

end