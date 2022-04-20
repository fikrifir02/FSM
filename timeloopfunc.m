%Function that runs time-loop for multiple species
%MPA is same size for all species and it is implemented at the same time

function xres = timeloopfunc(no_reefs, reefs_status,reefs_conn, reefs_std, umat, index_sp, Ninit,spparams, tmax, startmpa,Frem)

if isa(Ninit,'double') == 1
N1 = arrayfun(@(index_sp) preallfunc(Ninit, spparams, umat, no_reefs, reefs_std),index_sp, 'UniformOutput',false);
else
N1 = Ninit; %Standardized Ninitial by reef coverage
end
Biomass = nan(tmax, no_reefs);% copy nspp times for biomass in nmpa, reserve and fmpa for each species
Hvals = nan(tmax, 1); % total harvest in nmpa, reserve and fmpa for each species
Catches = nan(tmax, no_reefs);% catch in nmpa, reserve and fmpa for each species

for t  =  2:tmax
        if t==startmpa 
        %implement MPA
	    %***Turn this on or off to affect effort aggregation
        % Calculate number of reefs allocated as reserve, FMPA, and NMPA
        reef_reserve_no = numel(find(reefs_status(t,:)==1)); % number of reefs with reserve status
        reef_fmpa_no = numel(find(reefs_status(t,:)==0.5)); %number of reefs with fmpa status
        reef_nmpa_no = no_reefs - reef_reserve_no - reef_fmpa_no; %number of reefs with NMPA status
        
        % effort in each individual reefs before MPA implementation
        effort_u = umat{1,1}(1,1); 
        
        %estimate effort/fishing effort displacement after MPA
        %establishment in MPA (assumed well distributed effort)
        
        displaced_u_fmpa = (1-Frem)*effort_u*reef_fmpa_no; % amount of pressure displacef from FMPA
        displaced_u_reserve = effort_u* reef_reserve_no; % amount of pressure displaced from reserve
        effort_u_nmpa = ((displaced_u_fmpa + displaced_u_reserve+ (effort_u*reef_nmpa_no)))/(reef_nmpa_no); % exist and reallocated fishing pressure in NMPA
        
            for i = 1: no_reefs
                if reefs_status (t,i) == 1
                umat{1,1}(:,i) = 0 *umat{1,1}(:,i); %fishing mortality in reserve
                elseif reefs_status(t,i) == 0.5
                umat{1,1}(:,i) = Frem*umat{1,1}(:,i); %fishing mortality in fmpa
                else
                umat{1,1}(:,i) = effort_u_nmpa; %fishing mortality in NMPA
                end
            end
        else
        end
  % Aging& growth
   N1{1,1}(2:spparams.nages,:) = N1{1,1}(1:(spparams.nages-1),:);
   N1{1,1}(1,:) = 0; %remove recruits
   N1{1,1};
 
harv = zeros(1);
catch_zones = zeros(1, no_reefs);
    for tmonths = 1:spparams.nmonths
	N1 = arrayfun(@(index_sp) popmod(spparams, umat,N1,spparams.dT,reefs_std),index_sp, 'UniformOutput',false); %abundance at month t
	[harv1] = (arrayfun(@(index_sp) fishmod(N1, umat,spparams, reefs_std ,spparams.dT),index_sp, 'UniformOutput',false)); %harvest at month t
    [catch_zones1] = (arrayfun(@(index_sp) fishmod_zones(N1, umat,spparams,spparams.dT),index_sp, 'UniformOutput',false)); %harvest at month t
    harv1 = cell2mat(squeeze(harv1));
    catch_zones1 = cell2mat(squeeze(catch_zones1));
    harv = harv + harv1; %total harvest in whole mont for each species, each harvest at month t is accumulated (harvest at year)
    catch_zones = catch_zones + catch_zones1;   
    end
    
[Biom] = arrayfun(@(index_sp) getbio(N1, spparams),index_sp, 'UniformOutput',false);
Biomass(t,:) = Biom{1,1}; %biomass at year t in each reefs
Catches(t,:) = catch_zones; %catch at year t in each reefs 
Hvals(t,:) = harv; %total harvest at year t
end
xres.N1 = N1; %N1 abundance of each ages at t year.
xres.Biomass = Biomass; %biomass in each zones  for each species
xres.Catches = Catches; %catch in each zones  for each species
xres.Hvals = Hvals; %total harvest from all zones for each species
xres.umatvals = umat(1,:);

end