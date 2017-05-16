function [rxns_idx,rxnNames] = findExchRxnsFromMets(model,mets)
%FINDEXCHRXNSFROMMETS Find exchange reactions from extracellular metabolite names
%
% [rxns_idx,rxn_names] = findExchRxnsFromMets(model,metNames)
%
%REQUIRED INPUTS
% The model structure must contain the following fields:
%   model.S:        Stoichiometric matrix
%   model.mets:     Metabolite names
%   model.rxnNames: Reaction names
% mets specifies the exchange metabolites of interest (cell string)
%
%OUTPUT
% rxns_idx is a numeric array with the indices of the reactions
% rxns is a cell array with the reaction names
%
% Meghan Thommes 05/16/2017 - Changed name
% Meghan Thommes 06/14/2016 - New method to find exchange metabolites
% Meghan Thommes 05/12/2015

% Metabolites
[~,~,mets_idx] = intersect(mets,model.mets,'stable'); % index in model.mets
rxns_idx = zeros(size(mets_idx));
for idx = 1:length(mets_idx)
    idx_neg1 = find(model.S(mets_idx(idx),:) == -1); % find where met coeff = -1
    % find where reaction only has 1 met
    for ii = 1:numel(idx_neg1)
        num_met(ii) = numel(find(model.S(:,idx_neg1(ii)) ~= 0));
    end
    rxns_idx(idx) = idx_neg1(num_met == 1); % find column with only coeff = -1
    clear idx_neg1 num_met
end

rxnNames = model.rxnNames(rxns_idx);

end



