function [trspt_rxns] = identifyTrsptRxns(model)
%%identifyTrsptRxns Identify transport reactions (extracellular to
%%intracellular only).
%   trspt_rxns = identifyTrsptRxns(model)
%
%REQUIRED INPUT
% The model structure must contain the following fields:
%   model.S: Stoichiometric matrix
%
%OUTPUT
% trspt_rxns is a vector indicating the indices of the transport reactions
%   in the model
%
% Meghan Thommes 07/14/2017

%% Check Inputs

if (nargin < 1)
    error('myfuns:identifyTrsptRxns:NotEnoughInputs', ...
        'Not enough inputs: need a model file');
else
    if ~isstruct(model)
        error('myfuns:identifyTrsptRxns:IncorrectInput', ...
            '"model" needs to be a structure');
    elseif ~isfield(model,'S')
        error('myfuns:identifyTrsptRxns:IncorrectInput', ...
            '"model" needs "S" field');
    end
end

%% Identify Transport Reactions

% Identify Exchange Reactions
[exch_rxns,~] = identifyExchRxns(model);
% Identify Extracellular Metabolites
[exch_mets,~] = identifyExchMets(model);

% Identify Transport Reactions
trspt_rxns = [];
for ii = 1:numel(exch_mets)
    [~,C_0] = find(model.S(exch_mets(ii),:));
    trspt_rxns = [trspt_rxns, setdiff(C_0,exch_rxns)];
end
trspt_rxns = sort(unique(trspt_rxns))';

if numel(trspt_rxns) < numel(exch_rxns)
    error('myfuns:identifyTrsptRxns:IncorrectCalc', ...
        'Number of trspt_rxns is less than the number of exch_rxns');
end

end