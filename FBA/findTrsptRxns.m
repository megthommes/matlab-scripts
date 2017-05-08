function [trspt_rxns] = findTrsptRxns(model)
%%FINDTRSPTRXNS Find exchange and "nutrient" reactions
%   trspt_rxns = FINDTRSPTRXNS(model)
%
%
%REQUIRED INPUT
% The model structure must contain the following fields:
%   model.S:     Stoichiometric matrix
%
%OUTPUT
% trspt_rxns is a vector indicating the indices of the transport reactions
%   in the model
%
% Meghan Thommes 03/02/2017

%%

% Check Inputs
if (nargin < 1)
    error('Not enough inputs: need a model file');
else
    if ~isstruct(model)
        error('"model" needs to be a structure');
    elseif ~isfield(model,'S')
        error('"model" needs "S" field');
    end
end

% Identify Exchange Reactions
[exch_rxns,~] = findExchRxns(model);
exch_rxns = find(exch_rxns ~= 0);
% Identify Extracellular Metabolites
[exch_mets,~] = findExchMets(model);

% Identify Transport Reactions
trspt_rxns = [];
for ii = 1:numel(exch_mets)
    idx = find(full(model.S(exch_mets(ii),:))~=0);
    idx = setdiff(idx,exch_rxns);
    trspt_rxns = [trspt_rxns, idx];
end
trspt_rxns = sort(unique(trspt_rxns))';

end