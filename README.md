# Com_model
Run this code to make your community model with Matlab

Dependencies: Matlab_R2018b, ILOG_Cplex v.12.8, at least one model in .mat format 
Specify their locations at the begining of the Two_Org_Ratio

To run the code run Two_org_ratio, which depends on the function Make_Ratio_Com
This code works with all BiGG format GSM and with the new yeast model from Charlmers https://github.com/SysBioChalmers/yeast-GEM

In Two_org_ratio specify if you want amino acids to be exchangeable, if you want equal growth rate or competition for resources using 0 or 1 
Alternatively, ratios of communities can be chosen, transforming the problem's units to mmol-metabolites/(g-Total biomass DW * h).

After creating the community model in BiGG notation, an FBA will be done
For any questions
alexandre.tremblay@mail.utoronto.ca

Alexandre
