function [rxns_idx] = findExchRxnsFromExchMets(model,mets_idx,revFlag)
%findExchRxnsFromExchMets Find exchange reactions from extracellular 
%metabolite indices
%   rxns_idx = findExchRxnsFromExchMets(model,mets_idx)
%   rxns_idx = findExchRxnsFromExchMets(model,mets_idx,revFlag)
%
%By default, findExchRxnsFromExchMets assumes that exchange reactions are
% written as export reactions: 'A <==>'; therefore, uptaking a metabolite is
% a negative flux and the stoichiometric coefficient of metabolite A is -1.
%   If this is not the case, set revFlag=1.
%
%REQUIRED INPUTS
% The model structure must contain the following field:
%   model.S:        Stoichiometric matrix
% mets_idx: Vector with the indices of the exchange metabolites
%
%OPTIONAL INPUT
% revFlag: Indicates if the exchange reactions are reversed ("1") or not ("0")
%
%OUTPUTS
% rxns_idx: Vector with the indices of the exchange reactions
%
% Meghan Thommes 07/14/2017

%% Check Inputs

if (nargin < 2)
    error('myfuns:findExchRxnsFromExchMets:NotEnoughInputs', ...
        'Not enough inputs');
elseif (nargin == 2)
    if ~isstruct(model)
        error('myfuns:findExchRxnsFromExchMets:IncorrectInput', ...
            '"model" needs to be a structure');
    elseif ~isfield(model,'S')
        error('myfuns:findExchRxnsFromExchMets:IncorrectInput', ...
            '"model" needs "S"');
    end
    revFlag = 0;
else
    if ~isstruct(model)
        error('myfuns:findExchRxnsFromExchMets:IncorrectInput', ...
            '"model" needs to be a structure');
    elseif ~isfield(model,'S')
        error('myfuns:findExchRxnsFromExchMets:IncorrectInput', ...
            '"model" needs "S" field');
    end
    if ~isscalar(revFlag)
        error('myfuns:findExchRxnsFromExchMets:IncorrectInput', ...
            '"revFlag" needs to be a scalar');
    elseif (revFlag ~= 0 && revFlag ~= 1)
        error('myfuns:findExchRxnsFromExchMets:IncorrectInput', ...
            '"revFlag" needs to be "0" or "1"');
    end
end

%% Find exchange reactions from extracellular metabolite indices

% Find Exchange Reactions
exch_rxns = findExchRxns(model);

if ~revFlag % exchange reactions are not reversed
    % Find Reactions (Columns) with Stoichiometry of -1
    [~,C_1] = find(model.S(mets_idx,:) == -1);
    C_1 = unique(C_1);
else % exchange reactions are reversed
    % Find Reactions (Columns) with Stoichiometry of 1
    [~,C_1] = find(model.S(mets_idx,:) ==  1);
    C_1 = unique(C_1);
end

% Find Exchange Reaction Indices Corresponding to Metabolites
rxns_idx = intersect(C_1,exch_rxns);

end



