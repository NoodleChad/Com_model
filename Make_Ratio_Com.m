function [Com,Org1,Org2]=Make_Ratio_Com(model1,model2,Ratio1,Ratio2,Equalbio,AA_Exchange,Competition)
%% Rename models to fit standard notation
%find if there is [e]
  [n,~]=size(model1.mets(contains(model1.mets,'[e]'))); 
%if there is then change it
  if n>0
    model1.mets=strrep(model1.mets,'-','__');
    model1.rxns=strrep(model1.rxns,'-','__');
    model1.mets=strrep(model1.mets,'[e]','_e');
    model1.mets=strrep(model1.mets,'[c]','_c');
    model1.rxns=strrep(model1.rxns,'(e)','_e');
    model1.rxns=strrep(model1.rxns,'EX_glc_e','EX_glc__D_e');
  end
%model 2
  [n,~]=size(model2.mets(contains(model2.mets,'[e]'))); 
%if there is then change it
  if n>0
    model2.mets=strrep(model2.mets,'-','__');
    model2.rxns=strrep(model2.rxns,'-','__');
    model2.mets=strrep(model2.mets,'[e]','_e');
    model2.mets=strrep(model2.mets,'[c]','_c');
    model2.rxns=strrep(model2.rxns,'(e)','_e');
    model2.rxns=strrep(model2.rxns,'EX_glc_e','EX_glc__D_e');
  end
% Sometimes, there is an s at MetCharge, sometimes not
if isfield(model1,'metCharges') == 1
    model1.metCharge=model1.metCharges;
end
if isfield(model2,'metCharges') == 1
    model2.metCharge=model2.metCharges;
end
% Sometimes, grRules is just rules
if isfield(model1,'rules') == 1
    model1.grRules=model1.rules;
end
if isfield(model2,'rules') == 1
    model2.grRules=model2.rules;
end
% Sometimes it does not exist but do not worry
if isfield(model1,'metCharge') == 0
    model1.metCharge=zeros(size(model1.mets,1),1);
end
if isfield(model2,'metCharge') == 0
    model2.metCharge=zeros(size(model2.mets,1),1);
end
%% Order models in alphabetical order

[model1.mets,sortMets]=sort(model1.mets);
model1.metCharge=model1.metCharge(sortMets);
model1.metNames=model1.metNames(sortMets);
model1.metFormulas=model1.metFormulas(sortMets);
model1.S=model1.S(sortMets,1:end);
[model1.rxns,sortRxns]=sort(model1.rxns);
model1.S=model1.S(1:end,sortRxns);
if isfield(model1,'rxnGeneMat') == 1
    if isfield(model2,'rxnGeneMat') == 1
        model1.rxnGeneMat=model1.rxnGeneMat(sortRxns,1:end);
    end
end
model1.grRules=model1.grRules(sortRxns);
model1.rxnNames=model1.rxnNames(sortRxns);
model1.lb=model1.lb(sortRxns);
model1.ub=model1.ub(sortRxns);
model1.c=model1.c(sortRxns);
if isfield(model1,'rev') == 1
    if isfield(model2,'rev') == 1
        model1.rev=model1.rev(sortRxns);
    end
end
[model2.mets,sortMets]=sort(model2.mets);
model2.metCharge=model2.metCharge(sortMets);
model2.metNames=model2.metNames(sortMets);
model2.metFormulas=model2.metFormulas(sortMets);
model2.S=model2.S(sortMets,1:end);
[model2.rxns,sortRxns]=sort(model2.rxns);
model2.S=model2.S(1:end,sortRxns);
if isfield(model1,'rxnGeneMat') == 1
    if isfield(model2,'rxnGeneMat') == 1
        model2.rxnGeneMat=model2.rxnGeneMat(sortRxns,1:end);
    end
end
model2.grRules=model2.grRules(sortRxns);
model2.rxnNames=model2.rxnNames(sortRxns);
model2.lb=model2.lb(sortRxns);
model2.ub=model2.ub(sortRxns);
model2.c=model2.c(sortRxns);
if isfield(model1,'rev') == 1
    if isfield(model2,'rev') == 1
        model2.rev=model2.rev(sortRxns);
    end
end

%% limit ub and lb to have the same
    Max = 100;
    model1.lb(model1.lb<-Max)=-Max;
    model1.ub(model1.ub>Max)=Max;
    model2.lb(model2.lb<-Max)=-Max;
    model2.ub(model2.ub>Max)=Max;

%% identify EX reaction and metabolites
%Find where are the exchanged metabolites and if their name is alright
for n1 = 1:size(model1.S,2)
    nbz(n1) = sum(model1.S(:,n1)~=0); 
    if nbz(n1) == 1 & isempty(find(contains(model1.rxnNames(n1),"exchange"))) ~= 1;
        EX(n1) = 1; % there will be a 1 if this rxns is an exchange
        emets(n1) = find(model1.S(:,n1)~=0); %location of ex_metabolites
    else
        EX(n1) = 0;
        emets(n1) = 0;
    end
end
emets_org1 = emets(find(~ismember(emets,0)));
EX_org1 = find(ismember(EX,1));
if isempty(find(~contains(model1.mets(emets_org1),'_e'))) ~=1
    eloc_org1 = find(~contains(model1.mets(emets_org1),'_e'));  
    extra_org1 = find(contains(model1.metNames(emets_org1(eloc_org1)),'extracellular'));
    model1.mets(emets_org1(eloc_org1(extra_org1))) = strcat(model1.mets(emets_org1(eloc_org1(extra_org1))),'_e');
end
model1.rxns(EX_org1(find(~contains(model1.rxns(EX_org1),'EX'))))=...
    strcat('EX_',model1.rxns(EX_org1(find(~contains(model1.rxns(EX_org1),'EX')))));

% Identify and change structure
    pattern_mets ="_e";
    Exit_mets = contains(model1.mets,pattern_mets);
    A=find(ismember(Exit_mets,1)); 
    pattern_rxns ="EX_";
    Exit_rxns=contains(model1.rxns,pattern_rxns);
    B=find(ismember(Exit_rxns,1));
    [EX_mets_m,~]=size(A);
    [EX_rxns_m,~]=size(B);
%% Manipulate the matrix to have external metabolites first
    [m,n]=size(model1.S);
    C=find(ismember(Exit_mets,0)); % get all the internal and transport
    D=find(ismember(Exit_rxns,0)); % get all the internal and transport
    S1=sparse(m,n);
    S1(1:EX_mets_m,1:n)=model1.S(A,1:n);
    S1(EX_mets_m+1:m,1:n)=model1.S(C,1:n);
    S2=sparse(m,n);
    S2(1:m,1:EX_rxns_m)=S1(1:m,B);
    S2(1:m,EX_rxns_m+1:n)=S1(1:m,D);
%% Manipulate the rest of the variables in the model to do the same
% Now change order
Org1.mets(1:EX_mets_m,1)=model1.mets(A);
Org1.mets(EX_mets_m+1:m,1)=model1.mets(C);
Org1.metNames(1:EX_mets_m,1)=model1.metNames(A);
Org1.metNames(EX_mets_m+1:m,1)=model1.metNames(C);
Org1.metFormulas(1:EX_mets_m,1)=model1.metFormulas(A);
Org1.metFormulas(EX_mets_m+1:m,1)=model1.metFormulas(C);
Org1.metCharge(1:EX_mets_m,1)=model1.metCharge(A);
Org1.metCharge(EX_mets_m+1:m,1)=model1.metCharge(C);
Org1.genes=model1.genes;
if isfield(model1,'rxnGeneMat') == 1
    if isfield(model2,'rxnGeneMat') == 1
        Org1.rxnGeneMat=model1.rxnGeneMat;
        Org1.rxnGeneMat(1:EX_rxns_m,1:end)=model1.rxnGeneMat(B,1:end);
        Org1.rxnGeneMat(EX_rxns_m+1:n,1:end)=model1.rxnGeneMat(D,1:end);
    end
end
Org1.grRules(1:EX_rxns_m,1)=model1.grRules(B);
Org1.grRules(EX_rxns_m+1:n,1)=model1.grRules(D);
Org1.rxns(1:EX_rxns_m,1)=model1.rxns(B);
Org1.rxns(EX_rxns_m+1:n,1)=model1.rxns(D);
Org1.rxnNames(1:EX_rxns_m,1)=model1.rxnNames(B);
Org1.rxnNames(EX_rxns_m+1:n,1)=model1.rxnNames(D);
Org1.S=S2;
clear S1
clear S2
Org1.lb(1:EX_rxns_m,1)=model1.lb(B);
Org1.lb(EX_rxns_m+1:n,1)=model1.lb(D);
Org1.ub(1:EX_rxns_m,1)=model1.ub(B);
Org1.ub(EX_rxns_m+1:n,1)=model1.ub(D);
Org1.c(1:EX_rxns_m,1)=model1.c(B);
Org1.c(EX_rxns_m+1:n,1)=model1.c(D);
Org1.b=zeros(m,1);
if isfield(model1,'rev') == 1
    if isfield(model2,'rev') == 1
        Org1.rev(1:EX_rxns_m,1)=model1.rev(B);
        Org1.rev(EX_rxns_m+1:n,1)=model1.rev(D);
    end
end
Org1.description=model1.description;


%% Second organism
clearvars EX emets n1 nbz
for n1 = 1:size(model2.S,2)
    nbz(n1) = sum(model2.S(:,n1)~=0); 
    if nbz(n1) == 1 & isempty(find(contains(model2.rxnNames(n1),"exchange"))) ~= 1;
        EX(n1) = 1; % there will be a 1 if this rxns is an exchange
        emets(n1) = find(model2.S(:,n1)~=0); %location of ex_metabolites
    else
        EX(n1) = 0;
        emets(n1) = 0;
    end
end
emets_org2 = emets(find(~ismember(emets,0)));
EX_org2 = find(ismember(EX,1));
% Add _e to metabolites
if isempty(find(~contains(model2.mets(emets_org2),'_e'))) ~=1
    eloc_org2 = find(~contains(model2.mets(emets_org2),'_e'));  
    extra_org2 = find(contains(model2.metNames(emets_org2(eloc_org2)),'extracellular'));
    model2.mets(emets_org2(eloc_org2(extra_org2))) = strcat(model2.mets(emets_org2(eloc_org2(extra_org2))),'_e');
end
% Add EX_ to rxns
model2.rxns(EX_org2(find(~contains(model2.rxns(EX_org2),'EX'))))=...
    strcat('EX_',model2.rxns(EX_org2(find(~contains(model2.rxns(EX_org2),'EX')))));
%% identify EX reaction and metabolites
pattern_mets ="_e";
Exit_mets = contains(model2.mets,pattern_mets);
A=find(ismember(Exit_mets,1)); 
pattern_rxns ="EX_";
Exit_rxns=contains(model2.rxns,pattern_rxns);
B=find(ismember(Exit_rxns,1));
[EX2_mets_m,~]=size(A);
[EX2_rxns_m,~]=size(B);

%% Manipulate the matrix to have external metabolites first
[m2,n2]=size(model2.S);
C=find(ismember(Exit_mets,0)); % get all the internal and transport
D=find(ismember(Exit_rxns,0)); % get all the internal and transport
S1=sparse(m2,n2);
S1(1:EX2_mets_m,1:n2)=model2.S(A,1:n2);
S1(EX2_mets_m+1:m2,1:n2)=model2.S(C,1:n2);
S2=sparse(m2,n2);
S2(1:m2,1:EX2_rxns_m)=S1(1:m2,B);
S2(1:m2,EX2_rxns_m+1:n2)=S1(1:m2,D);

%% Manipulate the rest of the variables in the model2 to do the same

Org2.mets(1:EX2_mets_m,1)=model2.mets(A);
Org2.mets(EX2_mets_m+1:m2,1)=model2.mets(C);
Org2.metNames(1:EX2_mets_m,1)=model2.metNames(A);
Org2.metNames(EX2_mets_m+1:m2,1)=model2.metNames(C);
Org2.metFormulas(1:EX2_mets_m,1)=model2.metFormulas(A);
Org2.metFormulas(EX2_mets_m+1:m2,1)=model2.metFormulas(C);
Org2.metCharge(1:EX2_mets_m,1)=model2.metCharge(A);
Org2.metCharge(EX2_mets_m+1:m2,1)=model2.metCharge(C);
Org2.genes=model2.genes;
if isfield(model1,'rxnGeneMat') == 1
    if isfield(model2,'rxnGeneMat') == 1
        Org2.rxnGeneMat=model2.rxnGeneMat;
        Org2.rxnGeneMat(1:EX2_rxns_m,1:end)=model2.rxnGeneMat(B,1:end);
        Org2.rxnGeneMat(EX2_rxns_m+1:n2,1:end)=model2.rxnGeneMat(D,1:end);
    end
end
Org2.grRules(1:EX2_rxns_m,1)=model2.grRules(B);
Org2.grRules(EX2_rxns_m+1:n2,1)=model2.grRules(D);
Org2.rxns(1:EX2_rxns_m,1)=model2.rxns(B);
Org2.rxns(EX2_rxns_m+1:n2,1)=model2.rxns(D);
Org2.rxnNames(1:EX2_rxns_m,1)=model2.rxnNames(B);
Org2.rxnNames(EX2_rxns_m+1:n2,1)=model2.rxnNames(D);
Org2.S=S2;
clear S1
clear S2
Org2.lb(1:EX2_rxns_m,1)=model2.lb(B);
Org2.lb(EX2_rxns_m+1:n2,1)=model2.lb(D);
Org2.ub(1:EX2_rxns_m,1)=model2.ub(B);
Org2.ub(EX2_rxns_m+1:n2,1)=model2.ub(D);
Org2.b=zeros(m2,1);
Org2.c(1:EX2_rxns_m,1)=model2.c(B);
Org2.c(EX2_rxns_m+1:n2,1)=model2.c(D);
if isfield(model1,'rev') == 1
    if isfield(model2,'rev') == 1
        Org2.rev(1:EX2_rxns_m,1)=model2.rev(B);
        Org2.rev(EX2_rxns_m+1:n2,1)=model2.rev(D);
    end
end
Org2.description=model2.description;
%% Mix both organisms together
%First we need to clear all repetitions in the organisms 
All_mets(1:EX_mets_m,1)=Org1.mets(1:EX_mets_m,1);
All_mets(EX_mets_m+1:EX_mets_m+EX2_mets_m,1)=Org2.mets(1:EX2_mets_m,1);
y=unique(All_mets);
for k=1:size(y,1)
  freq(k)=sum(ismember(All_mets,y(k,:)));
end
Rep_mets=ismember(freq,2);
Rep_mets_name=y(Rep_mets);
Nonrep_mets=ismember(freq,1);
Nonrep_mets_name=y(Nonrep_mets);
Locate_mets=find(ismember(All_mets,Rep_mets_name)); 

All_rxns(1:EX_rxns_m,1)=Org1.rxns(1:EX_rxns_m,1);
All_rxns(EX_rxns_m+1:EX_rxns_m+EX2_rxns_m,1)=Org2.rxns(1:EX2_rxns_m,1);
y=unique(All_rxns);
for k=1:size(y,1)
  frec(k)=sum(ismember(All_rxns,y(k,:)));
end
Rep_rxns=ismember(frec,2);
Rep_rxns_name=y(Rep_rxns);
Nonrep_rxns=ismember(frec,1);
Nonrep_rxns_name=y(Nonrep_rxns);
Locate_rxns=find(ismember(All_rxns,Rep_rxns_name)); 

%number of shared and unique external metabolites & location in the matrix
[o,~]=size(Locate_mets);
nEx_mets_shared=o/2;
[o,~]=size(Locate_rxns);
nEx_rxns_shared=o/2;
clear o
[All_mets_size,~]=size(All_mets); 
[All_rxns_size,~]=size(All_rxns);

Org1_locate_mets=find(ismember(Org1.mets,Rep_mets_name));
Org2_locate_mets=ismember(Org2.mets,Rep_mets_name);
Org1_locate_rxns=find(ismember(Org1.rxns,Rep_rxns_name));
Org2_locate_rxns=ismember(Org2.rxns,Rep_rxns_name);

Org1_other_mets=find(ismember(Org1.mets,Nonrep_mets_name));
[~,~]=size(Org1_other_mets);
Org2_other_mets=find(ismember(Org2.mets,Nonrep_mets_name));
[~,~]=size(Org2_other_mets);

Org1_other_rxns=find(ismember(Org1.rxns,Nonrep_rxns_name));
[~,~]=size(Org1_other_rxns);
Org2_other_rxns=find(ismember(Org2.rxns,Nonrep_rxns_name));
[nEx_rxns_unique2,~]=size(Org2_other_rxns);


%% Making the mixed matrix

%Join all organisms for the external metabolites
S1=sparse(m+m2-nEx_mets_shared,n+n2);
S1(1:nEx_mets_shared,1:n)=Org1.S(Org1_locate_mets,1:end); %Org 1 share mets
S1(1:nEx_mets_shared,n+1:n+n2)=Org2.S(Org2_locate_mets,1:end);
S1(nEx_mets_shared+1:EX_mets_m,1:n)=Org1.S(Org1_other_mets,1:end);
S1(EX_mets_m+1:All_mets_size-nEx_mets_shared,n+1:n+n2)=Org2.S(Org2_other_mets,1:end);

%Join the rest of the metabolites
S1(All_mets_size-nEx_mets_shared+1:All_mets_size-nEx_mets_shared+m-EX_mets_m,1:n)=Org1.S(EX_mets_m+1:end,1:end);
S1(All_mets_size-nEx_mets_shared+m-EX_mets_m+1:end,n+1:n+n2)=Org2.S(EX2_mets_m+1:end,1:end);
[S1_m,~]=size(S1);

%Delete all shared EX_reactions and put all EX_reactions next to each other
S2=sparse(S1_m,n+n2-nEx_rxns_shared);
S2(1:S1_m,1:nEx_rxns_shared)=S1(1:end, Org1_locate_rxns);
S2(1:S1_m,nEx_rxns_shared+1:EX_rxns_m)=S1(1:end, Org1_other_rxns);
S2(1:S1_m,EX_rxns_m+1:All_rxns_size-nEx_rxns_shared)=S1(1:end, n+Org2_other_rxns);

S2(1:S1_m,All_rxns_size-nEx_rxns_shared+1:n+nEx_rxns_unique2)=S1(1:end,EX_rxns_m+1:n);
S2(1:S1_m,n+nEx_rxns_unique2+1:end)=S1(1:end,n+1+EX2_rxns_m:end);

%Start building the community model
Com.mets=Org1.mets(Org1_locate_mets);
Com.mets(nEx_mets_shared+1:EX_mets_m)=Org1.mets(Org1_other_mets);
Com.mets(EX_mets_m+1:All_mets_size-nEx_mets_shared)=Org2.mets(Org2_other_mets);
Com.mets(All_mets_size-nEx_mets_shared+1:All_mets_size-nEx_mets_shared+m-EX_mets_m)=Org1.mets(EX_mets_m+1:end);
Com.mets(All_mets_size-nEx_mets_shared+m-EX_mets_m+1:m+m2-nEx_mets_shared)=Org2.mets(EX2_mets_m+1:end);

Com.metNames=Org1.metNames(Org1_locate_mets);
Com.metNames(nEx_mets_shared+1:EX_mets_m)=Org1.metNames(Org1_other_mets);
Com.metNames(EX_mets_m+1:All_mets_size-nEx_mets_shared)=Org2.metNames(Org2_other_mets);
Com.metNames(All_mets_size-nEx_mets_shared+1:All_mets_size-nEx_mets_shared+m-EX_mets_m)=Org1.metNames(EX_mets_m+1:end);
Com.metNames(All_mets_size-nEx_mets_shared+m-EX_mets_m+1:m+m2-nEx_mets_shared)=Org2.metNames(EX2_mets_m+1:end);

Com.metFormulas=Org1.metFormulas(Org1_locate_mets);
Com.metFormulas(nEx_mets_shared+1:EX_mets_m)=Org1.metFormulas(Org1_other_mets);
Com.metFormulas(EX_mets_m+1:All_mets_size-nEx_mets_shared)=Org2.metFormulas(Org2_other_mets);
Com.metFormulas(All_mets_size-nEx_mets_shared+1:All_mets_size-nEx_mets_shared+m-EX_mets_m)=Org1.metFormulas(EX_mets_m+1:end);
Com.metFormulas(All_mets_size-nEx_mets_shared+m-EX_mets_m+1:m+m2-nEx_mets_shared)=Org2.metFormulas(EX2_mets_m+1:end);

Com.metCharge=Org1.metCharge(Org1_locate_mets);
Com.metCharge(nEx_mets_shared+1:EX_mets_m)=Org1.metCharge(Org1_other_mets);
Com.metCharge(EX_mets_m+1:All_mets_size-nEx_mets_shared)=Org2.metCharge(Org2_other_mets);
Com.metCharge(All_mets_size-nEx_mets_shared+1:All_mets_size-nEx_mets_shared+m-EX_mets_m)=Org1.metCharge(EX_mets_m+1:end);
Com.metCharge(All_mets_size-nEx_mets_shared+m-EX_mets_m+1:m+m2-nEx_mets_shared)=Org2.metCharge(EX2_mets_m+1:end);

[gene_org1,~]=size(Org1.genes);
[gene_org2,~]=size(Org2.genes);

Com.genes=Org1.genes;
Com.genes(gene_org1+1:gene_org1+gene_org2)=Org2.genes;

if isfield(model1,'rxnGeneMat') == 1
    if isfield(model2,'rxnGeneMat') == 1
        Com.rxnGeneMat=Org1.rxnGeneMat(Org1_locate_rxns,1:end);
        Com.rxnGeneMat(1:nEx_rxns_shared,gene_org1+1:gene_org2+gene_org1)=Org2.rxnGeneMat(Org2_locate_rxns,1:end);
        Com.rxnGeneMat(nEx_rxns_shared+1:EX_rxns_m,1:gene_org1)=Org1.rxnGeneMat(Org1_other_rxns,1:end);
        Com.rxnGeneMat(EX_rxns_m+1:All_rxns_size-nEx_rxns_shared,gene_org1+1:gene_org2+gene_org1)=Org2.rxnGeneMat(Org2_other_rxns,1:end);
        Com.rxnGeneMat(All_rxns_size-nEx_rxns_shared+1:n+nEx_rxns_unique2,1:gene_org1)=Org1.rxnGeneMat(EX_rxns_m+1:end,1:end);
        Com.rxnGeneMat(n+nEx_rxns_unique2+1:n+n2-nEx_rxns_shared,gene_org1+1:gene_org2+gene_org1)=Org2.rxnGeneMat(EX2_rxns_m+1:end,1:end);
    end
end
Com.grRules=Org1.grRules(Org1_locate_rxns);
Com.grRules(nEx_rxns_shared+1:EX_rxns_m)=Org1.grRules(Org1_other_rxns);
Com.grRules(EX_rxns_m+1:All_rxns_size-nEx_rxns_shared)=Org2.grRules(Org2_other_rxns);
Com.grRules(All_rxns_size-nEx_rxns_shared+1:n+nEx_rxns_unique2)=Org1.grRules(EX_rxns_m+1:end);
Com.grRules(n+nEx_rxns_unique2+1:n+n2-nEx_rxns_shared)=Org2.grRules(EX2_rxns_m+1:end);

Com.rxns=Org1.rxns(Org1_locate_rxns);
Com.rxns(nEx_rxns_shared+1:EX_rxns_m)=Org1.rxns(Org1_other_rxns);
Com.rxns(EX_rxns_m+1:All_rxns_size-nEx_rxns_shared)=Org2.rxns(Org2_other_rxns);
Com.rxns(All_rxns_size-nEx_rxns_shared+1:n+nEx_rxns_unique2)=Org1.rxns(EX_rxns_m+1:end);
Com.rxns(n+nEx_rxns_unique2+1:n+n2-nEx_rxns_shared)=Org2.rxns(EX2_rxns_m+1:end);

Com.rxnNames=Org1.rxnNames(Org1_locate_rxns);
Com.rxnNames(nEx_rxns_shared+1:EX_rxns_m)=Org1.rxnNames(Org1_other_rxns);
Com.rxnNames(EX_rxns_m+1:All_rxns_size-nEx_rxns_shared)=Org2.rxnNames(Org2_other_rxns);
Com.rxnNames(All_rxns_size-nEx_rxns_shared+1:n+nEx_rxns_unique2)=Org1.rxnNames(EX_rxns_m+1:end);
Com.rxnNames(n+nEx_rxns_unique2+1:n+n2-nEx_rxns_shared)=Org2.rxnNames(EX2_rxns_m+1:end);

Com.S=S2;

Com.ub=Org1.ub(Org1_locate_rxns);
Com.ub(nEx_rxns_shared+1:EX_rxns_m)=Org1.ub(Org1_other_rxns);
Com.ub(EX_rxns_m+1:All_rxns_size-nEx_rxns_shared)=Org2.ub(Org2_other_rxns);
Com.ub(All_rxns_size-nEx_rxns_shared+1:n+nEx_rxns_unique2)=Org1.ub(EX_rxns_m+1:end);
Com.ub(n+nEx_rxns_unique2+1:n+n2-nEx_rxns_shared)=Org2.ub(EX2_rxns_m+1:end);

Com.lb=Org1.lb(Org1_locate_rxns);
Com.lb(nEx_rxns_shared+1:EX_rxns_m)=Org1.lb(Org1_other_rxns);
Com.lb(EX_rxns_m+1:All_rxns_size-nEx_rxns_shared)=Org2.lb(Org2_other_rxns);
Com.lb(All_rxns_size-nEx_rxns_shared+1:n+nEx_rxns_unique2)=Org1.lb(EX_rxns_m+1:end);
Com.lb(n+nEx_rxns_unique2+1:n+n2-nEx_rxns_shared)=Org2.lb(EX2_rxns_m+1:end);

[x,~]=size(Com.S);
Com.b=zeros(x,1);

Com.c=Org1.c(Org1_locate_rxns);
Com.c(nEx_rxns_shared+1:EX_rxns_m)=Org1.c(Org1_other_rxns);
Com.c(EX_rxns_m+1:All_rxns_size-nEx_rxns_shared)=Org2.c(Org2_other_rxns);
Com.c(All_rxns_size-nEx_rxns_shared+1:n+nEx_rxns_unique2)=Org1.c(EX_rxns_m+1:end);
Com.c(n+nEx_rxns_unique2+1:n+n2-nEx_rxns_shared)=Org2.c(EX2_rxns_m+1:end);

if isfield(model1,'rev') == 1
    if isfield(model2,'rev') == 1
        Com.rev=Org1.rev(Org1_locate_rxns);
        Com.rev(nEx_rxns_shared+1:EX_rxns_m)=Org1.rev(Org1_other_rxns);
        Com.rev(EX_rxns_m+1:All_rxns_size-nEx_rxns_shared)=Org2.rev(Org2_other_rxns);
        Com.rev(All_rxns_size-nEx_rxns_shared+1:n+nEx_rxns_unique2)=Org1.rev(EX_rxns_m+1:end);
        Com.rev(n+nEx_rxns_unique2+1:n+n2-nEx_rxns_shared)=Org2.rev(EX2_rxns_m+1:end);
    end
end

Com.description='Community';
save ('FVA_model.mat','Com')
save ('Org1.mat','Org1')
save ('Org2.mat','Org2')

%% Find where are the uptake from the environment to the organism
    % If we let every transport open, the organism will be able to ferment
    % eveything that can make pyruvate (Ethanol, acetate, adenine, etc.)
    % Thus, we will close those reactions.
    
clearvars -except Ratio1 Ratio2 Equalbio model1 model2 AA_Exchange Competition
load 'FVA_model.mat'
load 'Org1'
load 'Org2'

% Exchange rxns are the only one that only contains a single value in their
% columns

for n1 = 1:size(Com.S,2)
    nbz(n1) = sum(Com.S(:,n1)~=0); 
    if nbz(n1) == 1;
        EX(n1) = 1; % there will be a 1 if this rxns is an exchange
        emets(n1) = find(Com.S(:,n1)~=0); %location of ex_metabolites
    else
        EX(n1) = 0;
        emets(n1) = 0;
    end
end
emets = emets(find(~ismember(emets,0)));
EX = find(ismember(EX,1));
% find all rxns that contains this the external metabolite
for i = 1:size(emets,2)
    trs(i,1:size(find(~ismember(Com.S(emets(i),:),0)),2))...
        = find(~ismember(Com.S(emets(i),:),0));
    trs(i,size(find(~ismember(Com.S(emets(i),:),0)),2)+1:10) = 0;
end
% Find the EX rxns associated with the trs, 
EXnTrs = [0 0];
for i = 1:size(trs,1)
    for i2 = 1:size(nonzeros(trs(i,:)),1)-1
        EXnTrs(end+1,1) = trs(i,find(ismember(trs(i,:),EX)));
        TrsinRow = nonzeros(trs(i,find(~ismember(trs(i,:),EX))));
        EXnTrs(end,2) = TrsinRow(i2); % [EX_loc, tr_loc]
    end
end
EXnTrs(1,:) = []; 

% Clean reactions which only contains extracellular metabolites
for i3 = 1:size(EXnTrs,1)
    MetinRxn(1,1:size(find(~ismember(Com.S(:,EXnTrs(i3,2)),0)),1)) =...
        find(~ismember(Com.S(:,EXnTrs(i3,2)),0))';
    Rem(i3,1) = isempty(find(~ismember(nonzeros(MetinRxn(1,:)),emets)));  
end
EXnTrs(find(ismember(Rem,1)),:) = [];

% Many exchange with water and H+
EXnTrs(find(ismember(EXnTrs(:,1),find(ismember(Com.rxns,{'EX_h_e','EX_h2o_e'})))),:) = [];

% find the sign in the matrix
for i4 = 1:size(EXnTrs,1)
    Sign(i4,1) = Com.S(find(~ismember(Com.S(:,EXnTrs(i4,1)),0)),EXnTrs(i4,1));
    Sign(i4,2) = Com.S(find(~ismember(Com.S(:,EXnTrs(i4,1)),0)),EXnTrs(i4,2));
    if Sign(i4,1) + Sign(i4,2) < 0
        if Com.ub(EXnTrs(i4,2)) ~= 0 
            if Com.lb(EXnTrs(i4,1))~=Com.ub(EXnTrs(i4,1))
                   Com.ub(EXnTrs(i4,2)) = -Com.lb(EXnTrs(i4,1))/abs(Sign(i4,2));
            end
        end
    elseif Sign(i4,1) + Sign(i4,2) >= 0
        if Com.lb(EXnTrs(i4,2)) ~= 0
            if Com.lb(EXnTrs(i4,1))~=Com.ub(EXnTrs(i4,1))
                   Com.lb(EXnTrs(i4,2)) = Com.lb(EXnTrs(i4,1))/abs(Sign(i4,2));
            end
        end
    else
       disp('Error with transport')
    end
end


Trs = Com.rxns(EXnTrs(:,2));
save ('transport.mat','Trs');
%% Now everything is closed, but sometimes it can be exchanged
% Such is the case of vitamins and amino acids
if AA_Exchange == 1
    Amino={'ALAtex','ARGtex','ASPtex','CYStex','GLNtex','GLUtex','GLYtex',...
        'HIStex','ILEtex','LEUtex','LYStex','METtex','PHEtex','PROtex',...
        'SERtex','THRtex','TRPtex','TYRtex','TYRtex','VALtex'};
    Amino_pos = find(ismember(Com.rxns,Amino'));
    Com.lb(Amino_pos)=-100;
    Com.ub(Amino_pos)=100;
    Vitamin={'THYMtex'};
    Vita_pos = find(ismember(Com.rxns,Vitamin));
    Com.lb(Vita_pos)=-100;
    Com.ub(Vita_pos)=100;
    if isfield(model1,'open') == 1
        Open_pos = find(ismember(Com.rxns,model1.open));
        Com.lb(Open_pos)=-100;
        Com.ub(Open_pos)=100;
    end
end

%% This step will devide the transport metabolites in each reaction
% This transformation stupulate that for one mmol glucose/(h*gDW) we have
% 1/Ratio inside the organism. This is done because gDW is different inside
% one organism and the whole community. 1 gDW com is 1*Ratio gDW of org 1

if Competition == 0
pattern_rxns ="EX_";
pattern_mets ="_e";
[Ex_rxns_Com,~]=size(Com.rxns(contains(Com.rxns,pattern_rxns)));
[Ex_rxns_Org1,~]=size(Com.rxns(contains(Org1.rxns,pattern_rxns)));
[Ex_rxns_Org2,~]=size(Com.rxns(contains(Org2.rxns,pattern_rxns)));
[emets_Com,~]=size(Com.mets(contains(Com.mets,pattern_mets)));
[emets_Org1,~]=size(Com.mets(contains(Org1.mets,pattern_mets)));
[emets_Org2,~]=size(Com.mets(contains(Org2.mets,pattern_mets)));

[Org1_rxns_size,~]=size(Org1.rxns);
[Org2_rxns_size,~]=size(Org2.rxns);
[Org1_mets_size,~]=size(Org1.mets);
[Org2_mets_size,~]=size(Org2.mets);
mets_size=Org1_mets_size+Org2_mets_size+emets_Com-emets_Org1-emets_Org2;
rxns_size=Org1_rxns_size+Org2_rxns_size+Ex_rxns_Com-Ex_rxns_Org1-Ex_rxns_Org2;
    
% Change the value of the ratio
    % For Org1
        Com.S(emets_Com+1:Org1_mets_size+(emets_Com-emets_Org1),Ex_rxns_Com+1:Ex_rxns_Org1+(Ex_rxns_Com-Ex_rxns_Org1))=Com.S(emets_Com+1:Org1_mets_size+(emets_Com-emets_Org1),Ex_rxns_Com+1:Ex_rxns_Org1+(Ex_rxns_Com-Ex_rxns_Org1))*(Ratio1);

    % For Org2
        Com.S(Org1_mets_size+(emets_Com-emets_Org1+1):mets_size,Org1_rxns_size+(Ex_rxns_Com-Ex_rxns_Org1)+1:rxns_size)=Com.S(Org1_mets_size+(emets_Com-emets_Org1+1):mets_size,Org1_rxns_size+(Ex_rxns_Com-Ex_rxns_Org1)+1:rxns_size)*(Ratio2);
else
    pattern_rxns ="EX_";
    pattern_mets ="_e";
    [Ex_rxns_Com,~]=size(Com.rxns(contains(Com.rxns,pattern_rxns)));
    [Ex_rxns_Org1,~]=size(Com.rxns(contains(Org1.rxns,pattern_rxns)));
    [Ex_rxns_Org2,~]=size(Com.rxns(contains(Org2.rxns,pattern_rxns)));
    [emets_Com,~]=size(Com.mets(contains(Com.mets,pattern_mets)));
    [emets_Org1,~]=size(Com.mets(contains(Org1.mets,pattern_mets)));
    [emets_Org2,~]=size(Com.mets(contains(Org2.mets,pattern_mets)));

    [Org1_rxns_size,~]=size(Org1.rxns);
    [Org2_rxns_size,~]=size(Org2.rxns);
    [Org1_mets_size,~]=size(Org1.mets);
    [Org2_mets_size,~]=size(Org2.mets);
    mets_size=Org1_mets_size+Org2_mets_size+emets_Com-emets_Org1-emets_Org2;
    rxns_size=Org1_rxns_size+Org2_rxns_size+Ex_rxns_Com-Ex_rxns_Org1-Ex_rxns_Org2;
end

%% Proportional biomass constraint 
% If both cells are growing at the same rate, then the steady state
% assumption is not violated

if Equalbio == 1
    obj = find(ismember(Com.c,1));
    Com.S(size(Com.S,1)+1,1:size(Com.S,2)) = 0; % Equal biomass constraint
    Com.S(size(Com.S,1),obj(1)) = 1/Ratio1;
    Com.S(size(Com.S,1),obj(2)) = -1/Ratio2;
    Com.mets(size(Com.S,1),1)={'Eq_biomass'};
    Com.metNames(size(Com.S,1),1)={'Eq_biomass'};
    Com.metFormulas(size(Com.S,1),1)={'Eq_biomass'};
    Com.metCharge(size(Com.S,1),1)=0;
    Com.b(size(Com.S,1),1)=0;
    
    %Changing the name on the model
    Org1_rxns = Com.rxns(Ex_rxns_Com+1:Org1_rxns_size-Ex_rxns_Org1+Ex_rxns_Com);
    Org2_rxns = Com.rxns(Org1_rxns_size-Ex_rxns_Org1+Ex_rxns_Com+1:rxns_size-1);
    Org1_rxns = strcat(Org1_rxns,'_org1');
    Org2_rxns = strcat(Org2_rxns,'_org2');
    Com.rxns(Ex_rxns_Com+1:Org1_rxns_size-Ex_rxns_Org1+Ex_rxns_Com)=Org1_rxns;
    Com.rxns(Org1_rxns_size-Ex_rxns_Org1+Ex_rxns_Com+1:size(Com.rxns,1)-1)=Org2_rxns;

elseif Equalbio == 0
    obj = find(ismember(Com.c,1));
    Com.S(size(Com.S,1)+1,1:size(Com.S,2)+1) = 0;
    Com.S(size(Com.S,1),size(Com.S,2)) = 1;
    Com.S(size(Com.S,1),obj(1)) = -1;
    Com.S(size(Com.S,1),obj(2)) = -1;
    Com.mets(size(Com.S,1),1)={'Eq_biomass'};
    Com.metNames(size(Com.S,1),1)={'Eq_biomass'};
    Com.metFormulas(size(Com.S,1),1)={'Eq_biomass'};
%     Com.metCharge(size(Com.S,1),1)=0;
    Com.b(size(Com.S,1),1)=0;
    if isfield(Com,'rev') == 1
        Com.rev(size(Com.S,2),1) = 0;
    end
    if isfield(Com,'rxnGeneMat') == 1
        Com.rxnGeneMat(1:size(Com.S,2),size(Com.rxnGeneMat,2)) = 0;
    end
    Com.rxnNames(size(Com.S,2),1)={'Eq_biomass'};
    Com.rxns(size(Com.S,2),1)={'Eq_biomass'};
    Com.grRules(size(Com.S,2),1)={''};
    Com.lb(size(Com.S,2),1) = 0;
    Com.ub(size(Com.S,2),1) = 100;
    Com.c(1:size(Com.S,2),1) = 0;
    Com.c(size(Com.S,2),1) = 1;
    Com.bio = [obj(1) obj(2)];
    %Changing the name on the model
    Org1_rxns = Com.rxns(Ex_rxns_Com+1:Org1_rxns_size-Ex_rxns_Org1+Ex_rxns_Com);
    Org2_rxns = Com.rxns(Org1_rxns_size-Ex_rxns_Org1+Ex_rxns_Com+1:rxns_size);
    Org1_rxns = strcat(Org1_rxns,'_org1');
    Org2_rxns = strcat(Org2_rxns,'_org2');
    Com.rxns(Ex_rxns_Com+1:Org1_rxns_size-Ex_rxns_Org1+Ex_rxns_Com)=Org1_rxns;
    Com.rxns(Org1_rxns_size-Ex_rxns_Org1+Ex_rxns_Com+1:size(Com.rxns,1)-1)=Org2_rxns;
end

%Add ratio in model
Com.Ratio=[Ratio1 Ratio2];

%% Ratio related LB and UB
% Lower bound should be affected by ratio too
if Competition == 0
   Com.ub(find(contains(Com.rxns,'_org1'))) =  Com.ub(find(contains(Com.rxns,'_org1'))) * Ratio1;
   Com.lb(find(contains(Com.rxns,'_org1'))) =  Com.lb(find(contains(Com.rxns,'_org1'))) * Ratio1;
   Com.ub(find(contains(Com.rxns,'_org2'))) =  Com.ub(find(contains(Com.rxns,'_org2'))) * Ratio2;
   Com.lb(find(contains(Com.rxns,'_org2'))) =  Com.lb(find(contains(Com.rxns,'_org2'))) * Ratio2;
end
%% Saving...
save ('Org1.mat','Org1')
save ('Org2.mat','Org2')
save ('2_organisms_model_ratio.mat','Com')
end