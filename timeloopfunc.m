%Function that runs time-loop for multiple species
%MPA is same size for all species and it is implemented at the same time

function xres = timeloopfunc(Aval_null, umat, index_sp, Ninit,spparams, Ainit, tmax, startmpa,Frem)

if isa(Ninit,'double') == 1
N1 = arrayfun(@(index_sp) preallfunc(Ninit, spparams, Ainit, umat),index_sp, 'UniformOutput',false);
else
N1 = Ninit;
end
Biomass = nan(tmax, 3);% copy nspp times for biomass in nmpa, reserve and fmpa for each species
Hvals = nan(tmax, 1); % total harvest in nmpa, reserve and fmpa for each species
Catches = nan(tmax, 3);% catch in nmpa, reserve and fmpa for each species

%browser()
A = Ainit; 

for t  =  2:tmax
        if t==startmpa && (Aval_null(1)+Aval_null(2))>0
        %implement MPA
         A(3) = Aval_null (2); %Aval_null(2) fmpa area proportion
         A(2) = Aval_null (1); %Aval_null(1) reserve area proportion
         A(1) = 1 - A(2) - A(3);
	    %***Turn this on or off to affect effort aggregation
		umat{1,1}(:,3) = Frem*umat{1,1}(:,3); %fishing mortality in fmpa
        umat{1,1}(:,1) = umat{1,1}(:,1) + umat{1,1}(:,1)*A(2) + (1-Frem)*umat{1,1}(:,1)*A(3); %fishing mortality in nmpa
        umat{1,1}(:,2) = 0;
        %***
          % half of the fishing pressure in the fmpa is allocated to nMPA
		%put fish into the MPA
            
		N1{1,1}(:,3) = N1{1,1}(:,1)*A(3);  %abundance in fmpa 
		N1{1,1}(:,2) = N1{1,1}(:,1)*A(2); %abundance in reserve 
        N1{1,1}(:,1) = N1{1,1}(:,1)*A(1); %abundance in nmpa 
        else
        end
  %Ageing/growth
   N1{1,1}(2:spparams.nages,:) = N1{1,1}(1:(spparams.nages-1),:);
   N1{1,1}(1,:) = 0; %remove recruits
   N1{1,1};
 
harv = zeros(1);
catch_zones = zeros(1, 3);
    for tmonths = 1:spparams.nmonths
	N1 = arrayfun(@(index_sp) popmod(spparams, umat, A, N1,spparams.dT),index_sp, 'UniformOutput',false); %abundance at month t
	[harv1] = (arrayfun(@(index_sp) fishmod(N1, umat,spparams, A ,spparams.dT),index_sp, 'UniformOutput',false)); %harvest at month t
    [catch_zones1] = (arrayfun(@(index_sp) fishmod_zones(N1, umat,spparams, A ,spparams.dT),index_sp, 'UniformOutput',false)); %harvest at month t
    harv1 = cell2mat(squeeze(harv1));
    catch_zones1 = cell2mat(squeeze(catch_zones1));
    harv = harv + harv1; %total harvest in whole mont for each species, each harvest at month t is accumulated (harvest at year)
    catch_zones = catch_zones + catch_zones1;   
    end
    
[Biom] = arrayfun(@(index_sp) getbio(N1, spparams),index_sp, 'UniformOutput',false);
Biomass(t,:) = cell2mat(Biom); %biomass at year t in each zones for each species
Catches(t,:) = catch_zones; %catch at year t in each zones for each species
Hvals(t,:) = harv; %total harvest at year t
end
xres.N1 = N1; %N1 abundance of each ages at t year.
xres.Biomass = Biomass; %biomass in each zones  for each species
xres.Catches = Catches; %catch in each zones  for each species
xres.Hvals = Hvals; %total harvest from all zones for each species
xres.umatvals = umat(1,:);

end