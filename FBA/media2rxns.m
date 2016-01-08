function [rxn_idx,rxn_name,met_idx,met_name] = media2rxns(model,media)
%media2rxns Find the exchange reactions associated with compounds
%
%REQUIRED INPUTS
% The model structure must contain the following fields:
%   model.S:     Stoichiometric matrix (m x n)
%   model.mets:  Metabolite IDs (m x 1)
%   model.rxns:  Reaction IDs (n x 1)
%
%OUTPUT
% rxn_idx:       Numeric array with the indices of reactions
% rxn_name:      Cell array with the reaction IDs
% met_idx:       Numeric array with the indices of metabolites
% met_name:      Cell array with the metabolite IDs
%
% 01/08/2016 Meghan Thommes

%% Check Inputs

%% Match Exchange Reactions to Extracellular Metabolites

% see if media metabolites are in model
modelMedia_idx = cellfun(@(x) strfind(model.mets,x),media, 'UniformOutput',false);

% media metabolites in the model
modelMedia = model.mets(cell2mat(cellfun(@(x) find(~cellfun('isempty',x)),modelMedia_idx, 'UniformOutput',false)));

% reactions that correspond to media metabolites
[rxn_idx_temp,~] = findRxnsFromMets(model,modelMedia);
rxn_idx = find(rxn_idx_temp == 1);
rxn_name = model.rxns(rxn_idx);

% media metabolites that have reactions
[met_idx,~] = findMetsFromRxns(model,rxn_idx);
met_name = model.mets(met_idx);

end