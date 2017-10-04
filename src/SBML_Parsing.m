
clear all
clc
model = 'BIOMD0000000050';
m1 = sbmlimport(strcat('C:/Users/sedaghat/Dropbox (Personal)/MONET_MetabolicNetworks/Data/RealMetabolicNetworks/', model,'.xml'));

[M,objSpecies,objReactions] = getstoichmatrix(m1);
Stoich = full(M);
reversibility = 1*cell2mat(get(get (m1, 'Reactions'),'Reversible'));

save(strcat('C:/Users/sedaghat/Dropbox (Personal)/MONET_MetabolicNetworks/Data/RealMetabolicNetworks/', model, '.mat'), 'Stoich', 'reversibility')

addpath('C:\Users\sedaghat\Dropbox (Personal)\MONET_MetabolicNetworks\Code\Other codes\efmtool');


