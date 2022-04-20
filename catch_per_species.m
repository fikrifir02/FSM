%function to accumulate catches for each species in all zones

function [catch_sp] = catch_per_species(xres)
            catch_sp = sum(xres.Catches(:,:),2);
end