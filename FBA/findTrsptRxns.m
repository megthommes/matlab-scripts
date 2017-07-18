function [trspt_rxns] = findTrsptRxns(model)
%%findTrsptRxns Find transport
%   trspt_rxns = findTrsptRxns(model)
%
%REQUIRED INPUT
% The model structure must contain the following fields:
%   model.S:     Stoichiometric matrix
%
%OUTPUT
% trspt_rxns is a vector indicating the indices of the transport reactions
%   in the model
%
% Meghan Thommes 07/14/2017
% Meghan Thommes 03/02/2017

%% Check Inputs
if (nargin < 1)
    error('Not enough inputs: need a model file');
else
    if ~isstruct(model)
        error('"model" needs to be a structure');
    elseif ~isfield(model,'S')
        error('"model" needs "S" field');
    end
end

%% Find Transport Reactions

% Identify Exchange Reactions
[exch_rxns,~] = findExchRxns(model);
% Identify Extracellular Metabolites
[exch_mets,~] = findExchMets(model);

% Identify Transport Reactions
trspt_rxns = [];
for ii = 1:numel(exch_mets)
    [~,C_0] = find(model.S(exch_mets(ii),:));
    trspt_rxns = [trspt_rxns, setdiff(C_0,exch_rxns)];
end
trspt_rxns = sort(unique(trspt_rxns))';

end