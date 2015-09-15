function [metInds,mets] = findMetsFromRxns(model,rxnInds)
%findMetsFromRxns Find metabolites from reaction indices - exchange
%reactions only
%
% [metInds,mets] = findMetsFromRxns(model,rxnInds)
%
%REQUIRED INPUTS
% The model structure must contain the following fields:
%   model.S:     Stoichiometric matrix
% rxnInds specifies the exchange reactions of interest. rxnInds can be
%   - logical array with the length of the reactions
%   - numeric array with the indices of the reactions
%
%OUTPUT
% metInds is a numeric array with the indices of the metabolites
% mets is a cell array with the metabolite names
%
%5/12/15 Meghan Thommes - Use reaction names instead of indices
%3/24/15 Meghan Thommes

exMets = model.S(:,rxnInds);
[metInds,~] = find(exMets ~= 0);
mets = model.metNames(metInds);

end

