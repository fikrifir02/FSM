% Non-spatial MPA model
% r (growth coeffient)is constant. Fishing skill varies
% Functions
% Brown 11 Nov 2013 - R version
% Fikri Nov 2019 - Matlab version
% V3 Has monthly time-steps
% V4 has multiple species
% V5 has adult movement
% Effort is aggregated


function [N] = popmod(spparams, umat, A, N1, dT)
    N = N1{1,1};
    Fmortmat = umat{1,1}*spparams.q*dT;
    %Movement, only if there is an MPA
     if (A(1)>0 && A(2)>0 && A(3)>0) 
     
     %Immigration and emigration, movement
     %nmpa (1), reserve (2), and fmpa (3)
     %from nmpa to reserve
     moverate_1_to_2 = A(2);
     %from nmpa to fmpa
     moverate_1_to_3 = A(3);
     %from reserve to fmpa
     moverate_2_to_3 = A(3);
     %from reserve to nmpa
     moverate_2_to_1 = A(1);
     %from fmpa to reserve
     moverate_3_to_2 = A(2);
     %from fmpa to nmpa
     moverate_3_to_1 = A(1);

   
     %movement
     move_1_to_2 = (spparams.mrate.* moverate_1_to_2*N(:,1));
     %from nmpa to fmpa
     move_1_to_3 = (spparams.mrate.* moverate_1_to_3*N(:,1));
     %from reserve to fmpa
     move_2_to_3 = (spparams.mrate.* moverate_2_to_3*N(:,2));
     %from reserve to nmpa
     move_2_to_1 = (spparams.mrate.* moverate_2_to_1*N(:,2));
     %from fmpa to reserve
     move_3_to_2 = (spparams.mrate.* moverate_3_to_2*N(:,3));
     %from fmpa to nmpa
     move_3_to_1 = (spparams.mrate.* moverate_3_to_1*N(:,3));
     
          
	 N(:,1) = N(:,1) - move_1_to_2 -  move_1_to_3 + move_2_to_1 + move_3_to_1; %number of fish remains in nmpa
	 N(:,2) = N(:,2) - move_2_to_3 -  move_2_to_1 + move_1_to_2 + move_3_to_2; %number of fish remains in reserve
     N(:,3) = N(:,3) - move_3_to_2 -  move_3_to_1 + move_1_to_3 + move_2_to_3; %number of fish remains in fmpa
    elseif A(1)== 0 %nmpa proportion equals to zero
     moverateN3 = A(2); 
     moverateN2 = A(3); 
     moverateN1 = 0; %zero because there is no proportion for nmpa
    
	 moveN1 = (spparams.mrate.* moverateN1*N(:,1)); %no movement from nmpa because there is no nmpa proportion
	 moveN2 = (spparams.mrate.* moverateN2*N(:,2)); 
     moveN3 = (spparams.mrate.* moverateN3*N(:,3)); 
     
	 N(:,1) = N(:,1) ; %number of fish remains in nmpa, just the original nmpa abundance which should be equal to zero
	 N(:,2) = N(:,2) - moveN2 + moveN3; 
     N(:,3) = N(:,3) + moveN2 - moveN3; 
     
     elseif A(2)== 0 %reserve proportion equals to zero
     moverateN3 = A(1); 
     moverateN2 = 0; %zero because there is no proportion for reserve
     moverateN1 = A(3);
    
	 moveN1 = (spparams.mrate.* moverateN1*N(:,1)); 
	 moveN2 = (spparams.mrate.* moverateN2*N(:,2)); %no movement from reserve because there is no reserve proportion
     moveN3 = (spparams.mrate.* moverateN3*N(:,3)); 
     
	 N(:,1) = N(:,1) - moveN1 + moveN3; 
	 N(:,2) = N(:,2) ; %number of fish remains in reserve, just the original reserve abundance which should  equal to zero
     N(:,3) = N(:,3) + moveN1 - moveN3; 
     
     elseif A(3)== 0 %fmpa proportion equals to zero
     moverateN3 = 0; %zero because there is no proportion for fmpa
     moverateN2 = A(1); 
     moverateN1 = A(2);
    
	 moveN1 = (spparams.mrate.* moverateN1*N(:,1)); 
	 moveN2 = (spparams.mrate.* moverateN2*N(:,2)); 
     moveN3 = (spparams.mrate.* moverateN3*N(:,3)); %no movement from fmpa because there is no reserve proportion
     
	 N(:,1) = N(:,1) - moveN1 + moveN2; 
	 N(:,2) = N(:,2) + moveN1 - moveN2;  
     N(:,3) = N(:,3) ;%number of fish remains in reserve, just the original reserve abundance which should  equal to zero
    else
    end
    %Mortality (from natural and fishing mortality)
    N(:,:) = N(:,:) .* exp((-(spparams.M*dT)) - Fmortmat(:,:));
   
    %Total larvae spawned
    spawn = sum(sum(N,2).* spparams.matind' .* spparams.weights'  .* spparams.fecparam'*dT);
    %Recruitment to each patch, relative to patch area
    relarea = [1 - A(2)-A(3), A(2),  A(3)];
	spawn_patch = spawn; % *.relarea %relarea is not used assuming spawn does not depend on spatial variable
	N(1,:) = N(1,:) + ((dT*spparams.alpha.* relarea .* spawn_patch) ./ (1 + spparams.beta * spawn_patch)); %put fish at age 1 (recruit), it is a density function on spawning biomass
	N(isnan(N)) = 0;
	N;
end

    