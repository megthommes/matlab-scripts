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
% Meghan Thommes 06/14/2016 - New method to find exchange metabolites
% Meghan Thommes 05/12/2015

% Metabolites
[~,~,mets_idx] = intersect(metNames,model.mets); % index in model.mets
rxns_idx = zeros(size(mets_idx));
for idx = 1:length(mets_idx)
    idx_neg1 = find(model.S(mets_idx(idx),:) == -1); % find where met coeff = -1
    rxns_idx(idx) = idx_neg1(sum(full(model.S(:,idx_neg1)),1) == -1); % find column with only coeff = -1
    clear idx_neg1
end

rxns = model.rxnNames(rxns_idx);

end



