% function [stoich_coeff,metID,metName] = rxnFormula(model,rxnID)

clear all
close all
clc

load('C:\Users\mthommes\OneDrive\Documents\BU\Lab\Metabolic_Models\Escherichia_coli\E_coli_Palsson.mat')
model = E_coli;
rxnID = model.rxns(7);

for rxn_num = 1:length(rxnID)
    [isRxn index] = ismember(rxnID{:},model.rxns); % find rxn index in model
    rxn_idx = index(isRxn);
end

% end