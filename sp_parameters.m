% Species parameters
% CJ Brown 9 Dec 2013 - R version
% parameters are modified to fit into three patch model
% F Firmansyah August 2020 - Matlab version

FPmult = 2;  %multiplier for how overfished it is (multiply by optimal effort)

nspp = num_sp; %number of species - FF
nmonths = 12; %time period, how annual growth is calclulated, divided by number of months - FF
dT = 1/nmonths;%time interval for calculating dyamic of abundance - FF
awght = 0.00001; %a-constant for calculating fish biomass - FF
bwght = 3; % b-constant for calculating fish biomas - FF

%minimum biomass (as prop of unfished) where profits are zero
minbio_prof = repelem(0.05, nspp); %minbio.prof changed to minbio_prof, otherwise object will be a struct object
%% Species 1 parrotfish Chlorurus sordidus
k = 1.083; %k - growth parameter - FF
linf = 193; %linf - l infinity - FF
maxage = 9; %maximum age
tzero = 0; %theoretical age when fish at zero cm length - FF
amat = 1; %age at maturity - FF
ages = 1:maxage; %fish age range - FF
awght1 = 1.82*(10^-5);
bwght1 = 3.15;

afish = 2; %age at vulnerable to fishery - FF
afishind = zeros(maxage, 3); %create nrow (maximum age), 3 columns matrix representing reserve, FMPA, nMPA - FF
afishind(ages>=amat,:) = 1; %recalculate afishind to determine fish vulnerable age at each zones - FF

lengths = linf *(1 - exp(-k * (ages - tzero))); %calculate length for each ages - FF
weights = awght1 * (lengths.^bwght1); %calculate weight for each ages - FF
weightsmat = [weights; weights; weights]'; %create an object which has three columns (nMPA, reserve, FMPA) and n rows (age range) -FF
matind = zeros(1, maxage); %prepare matrix for fish maturity 
matind(:,ages>=amat) = 1; %recalculate matind - fish older than or equals to amat = 1

%% Species 2 snapper: Lutjanus vittus
k2 = 0.37; 
linf2 = 325;
maxage2 = 12;
tzero2 = -0.23;
amat2 = 1;
ages2 = 1:maxage2;
awght2 = 9.99*(10^-6);
bwght2 = 3.086;

afish2 = 1;
afishind2 =zeros(maxage2, 3); 
afishind2(ages2>=amat2,:) = 1;

lengths2 = linf2 *(1 - exp(-k2 * (ages2 - tzero2)));
weights2 = awght2 * (lengths2.^bwght2);
weightsmat2 = [weights2; weights2; weights2]';  
matind2 = zeros(1, maxage2);
matind2(:,ages2>=amat2) = 1;

%% Species 3: Coral trout Plectropomus leopardus - see table 1 for parameter
k3 = 0.45;
linf3 = 522;
maxage3 = 18;
tzero3 = 0.77;
amat3 = 3;
ages3 = 1:maxage3;
awght3 = 7.8*(10^-3); %this param is in cm/kg, so need to divide length by 10
bwght3 = 3.157;

afish3 = 2;
afishind3 =zeros(maxage3, 3); 
afishind3(ages3>=amat3,:) = 1;

lengths3 = linf3 *(1 - exp(-k3 * (ages3 - tzero3)));
weights3 = awght3 * ((lengths3./10).^bwght3);
weightsmat3 = [weights3; weights3; weights3]';
matind3 = zeros(1, maxage3);
matind3(:,ages3>=amat3) = 1;

%% Species 4: Rabbit fish siganus argenteus - see table 1 for parameter
% based on https://www.publish.csiro.au/mf/Fulltext/MF16169 

k4 = 0.56; %k - growth parameter - FF
linf4 = 274; %linf - l infinity - FF
maxage4 = 8; %maximum age, rounded from 7.8
tzero4 = -0.30; %theoretical age when fish at zero cm length - FF
amat4 = 1.3; %age at maturity -for females
ages4 = 1:maxage4; %fish age range - FF
awght4 = 0.01440;
bwght4 = 3.238;

afish4 = 1; %age at vulnerable to fishery - FF
afishind4 = zeros(maxage4, 3); %create nrow (maximum age), 3 columns matrix representing reserve, FMPA, nMPA - FF
afishind4(ages4>=amat4,:) = 1; %recalculate afishind to determine fish vulnerable age at each zones - FF

lengths4 = linf3 *(1 - exp(-k3 * (ages4 - tzero4)));
weights4 = awght4 * ((lengths4./10).^bwght4);
weightsmat4 = [weights4; weights4; weights4]';
matind4 = zeros(1, maxage4);
matind4(:,ages4>=amat4) = 1;

%% Species parameterss
%compiling all parameters set above in a struct object
spparams = struct(); 
spparams.nages = [maxage, maxage2, maxage3, maxage4];
spparams.q = [1,1,1,1]; %catch coefficient
spparams.M = [0.31,0.34, 0.45, 0.56]; %natural mortality rate
spparams.matind = struct (); %matrix of individuals age
spparams.matind = {matind, matind2, matind3,matind4};
spparams.weights = struct ();
spparams.weights = {weights, weights2, weights3,weights4};
spparams.weightsmat = struct ();
spparams.weightsmat = {weightsmat, weightsmat2, weightsmat3, weightsmat4};
spparams.fecparam = [1, 1, 1, 1];
spparams.CR = [4, 4, 4, 4]; %what is CR? - FF
spparams.Rmax = [4, 1, 0.25, 16]; %maximum recruitment, used 4times magnitude, from low to high, siganus was the highest
spparams.afishind = struct ();
spparams.afishind = {afishind, afishind2, afishind3, afishind4};
spparams.RelVal = [1,2,4,1]; %wRelative value, value of rabbit fish and parrot fish are the same
spparams.alpha = [];% recruitment parameter
spparams.beta = []; % recruitment parameter
spparams.Lzero = []; % recruitment parameter
spparams.mrate = [0.01, 0.01, 0.01, 0.01]; %movement rate
spparams.omega = zeros(1,nspp); % cost parameter
spparams.minbio_prof = minbio_prof; %cost/profit parameter
spparams.nmonths = nmonths;
spparams.dT = dT;

for ispp = 1:nspp
	abund = spparams.Rmax(ispp).*exp(-spparams.M(ispp) .* (0:(spparams.nages(ispp)-1)));
	spparams.Lzero(ispp) = sum(abund .* spparams.matind{1,ispp} .* spparams.weights{1,ispp} .* spparams.fecparam(ispp).*dT);
	spparams.alpha(ispp) = (spparams.CR(ispp).*spparams.Rmax(ispp))/	spparams.Lzero(ispp);
	spparams.beta(ispp) = (spparams.CR(ispp)-1)/spparams.Lzero(ispp);
end