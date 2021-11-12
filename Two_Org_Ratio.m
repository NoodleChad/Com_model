clear;
cd /nfs/homes/tremb133/CommunityStuff/Ratio/model
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio');
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio/model');
% Cplex Path
addpath('/nfs/homes/tremb133/CPLEX_Studio128/cplex/matlab/x86-64_linux') 
%% Making the two models 

% Enter the model of your first organism here and rename
    load iML1515;
    model1=iML1515;
    Ratio1=0.5;
    
% Enter the model of your second organism here and rename
    load yeastGEM_v8_5;
    model2=model;
    Ratio2=1-Ratio1;
    
Equalbio = 0; %Equal biomass (1), unconstraned biomass (0)
AA_Exchange = 0; % Allow exchange of amino acids and vitamin between strains
Competition = 1; % There won't be any ratio and strain will compete for o2 and glc (1)

[Com,Org1,Org2]=...
    Make_Ratio_Com(model1,model2,Ratio1,Ratio2,Equalbio, AA_Exchange, Competition); % function

% Com.lb(find(contains(Com.rxns,'ATPM'))) = 0;
FBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
Biomass = FBA(find(ismember(Com.c,1)))
