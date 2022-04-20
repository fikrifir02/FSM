% Non-spatial MPA model
% r (growth coeffient)is constant. Fishing skill varies
% Functions
% Brown 11 Nov 2013 - R version
% Fikri Nov 2019 - Matlab version
% V3 Has monthly time-steps
% V4 has multiple species
% V5 has adult movement
% Effort is aggregated


function [N] = popmod(spparams, umat,N1, dT, reefs_std)
    N = N1{1,1}; 
    Fmortmat = umat{1,1}*spparams.q*dT;
    %Movement, only if there is an MPA
    
    
    % Fish movement
%      for i = 1: no_reefs
%          for so = 1: numel(reefs_conn(1,:))
%              for si = 1 : numel(reefs_conn(:,1))
%                  N(:,i) = N(:,i) + N(:,si).*reefs_conn(si,i).* spparams.mrate;    %reef i being sink of fish
%              end
%              N(:,i) = N(:,i) - N(:,i).*reefs_conn(i,so).* spparams.mrate;          %reef i being sink of fish
%          end
%      end
%             
    %Mortality (from natural and fishing mortality)
    N(:,:) = N(:,:) .* exp((-(spparams.M*dT)) - Fmortmat(:,:));
   
    %Total larvae spawned
    spawn = sum(sum(N,2).* spparams.matind' .* spparams.weights'  .* spparams.fecparam'*dT);
    %Recruitment to each patch, relative to patch area
	spawn_patch = spawn; % *.relarea %relarea is not used assuming spawn does not depend on spatial variable
	N(1,:) = N(1,:) + ((dT*spparams.alpha .* reefs_std'.*spawn_patch) ./ (1 + spparams.beta * spawn_patch)); %put fish at age 1 (recruit), it is a density function on spawning biomass
	N(isnan(N)) = 0;
	N;
end

    