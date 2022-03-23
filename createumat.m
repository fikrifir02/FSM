%Create fishing mortality list, umat, 
function[umat] = createumat(ispp,afishind,Btot) % this function is used to apply fishing mortality rate on each vulnerable ages
umat = Btot*afishind; %Btot represents fishing mortality rate, it is used on fishmod and popmod functions
end