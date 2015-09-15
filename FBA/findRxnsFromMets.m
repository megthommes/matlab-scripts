function [rxns_idx,rxns] = findRxnsFromMets(model,metNames)
%findRxnsFromMets Find reactions from metabolite names - exchange
%reactions only
%
% [rxnInds,rxns] = findRxnsFromMets(model,metNames)
%
%REQUIRED INPUTS
% The model structure must contain the following fields:
%   model.S:        Stoichiometric matrix
%   model.mets:     Metabolite names
%   model.rxnNames: Reaction names
% metNames specifies the exchange metabolites of interest (cell string)
%
%OUTPUT
% rxns_idx is a numeric array with the indices of the reactions
% rxns is a cell array with the reaction names
%
%5/12/15 Meghan Thommes

% Exchange Reactions
[excRxns,~] = findExchRxns(model);

% Metabolites
[isMet, index] = ismember(metNames,model.mets);
mets_idx = index(isMet)'; % index in model.mets
rxns_idx = full((sum(model.S(mets_idx,:)==-1,1) == 1) & (sum(model.S(mets_idx,:)~=0) == 1))';
rxns_idx(excRxns-rxns_idx ~= 0) = 0; % index in model.rxns

rxns = model.rxnNames(rxns_idx);

end

