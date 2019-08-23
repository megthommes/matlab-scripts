function [trspt_rxns] = findTrsptRxns(model,revFlag)
%%findTrsptRxns Find transport
%   trspt_rxns = findTrsptRxns(model)
%
%REQUIRED INPUT
% The model structure must contain the following fields:
%   model.S:     Stoichiometric matrix
%OPTIONAL INPUT
% revFlag: Indicates if the exchange reactions are reversed ("1") or not ("0")
%
%OUTPUT
% trspt_rxns is a vector indicating the indices of the transport reactions
%   in the model
%
% Meghan Thommes 07/14/2017
% Meghan Thommes 03/02/2017

%% Check Inputs
if (nargin < 1)
    error('myfuns:findTrsptRxns:NotEnoughInputs', ...
        'Not enough inputs: need a model file');
elseif (nargin == 1)
    if ~isstruct(model)
        error('myfuns:findTrsptRxns:NotEnoughInputs', ...
            '"model" needs to be a structure');
    elseif ~isfield(model,'S')
        error('myfuns:findTrsptRxns:NotEnoughInputs', ...
            '"model" needs "S" field');
    end
    revFlag = 0;
else
    if ~isstruct(model)
        error('myfuns:findTrsptRxns:NotEnoughInputs', ...
            '"model" needs to be a structure');
    elseif ~isfield(model,'S')
        error('myfuns:findTrsptRxns:NotEnoughInputs', ...
            '"model" needs "S" field');
    end
end

if (nargin > 1)
    if ~isscalar(revFlag)
        error('myfuns:findTrsptRxns:IncorrectInput', ...
            '"revFlag" needs to be a scalar');
    elseif (revFlag ~= 0 && revFlag ~= 1)
        error('myfuns:findTrsptRxns:IncorrectInput', ...
            '"revFlag" needs to be "0" or "1"');
    end
end

%% Find Transport Reactions

% Identify Exchange Reactions
[exch_rxns,~] = findExchRxns(model,revFlag);
% Identify Extracellular Metabolites
[exch_mets,~] = findExchMets(model,revFlag);

% Identify Transport Reactions
trspt_rxns = [];
for ii = 1:numel(exch_mets)
    [~,C_0] = find(model.S(exch_mets(ii),:));
    trspt_rxns = [trspt_rxns, setdiff(C_0,exch_rxns)];
end
trspt_rxns = sort(unique(trspt_rxns))';

end