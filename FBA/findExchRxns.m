function [exch_rxns,med_rxns] = findExchRxns(model,revFlag)
%%findExchRxns Find exchange reactions and identify "medium" reactions
%   [exch_rxns,med_rxns] = findExchRxns(model)
%   [exch_rxns,med_rxns] = findExchRxns(model,revFlag)
%
%By default, findExchRxns assumes that exchange reactions are written as
%export reactions: 'A <==>'; therefore, uptaking a metabolite is a
%negative flux and the stoichiometric coefficient of metabolite A is -1.
%   If this is not the case, set revFlag=1.
%
%REQUIRED INPUT
% The model structure must contain the following fields:
%   model.S:     Stoichiometric matrix
%
%OPTIONAL INPUTS
% The model structure my contain the following field:
%   model.c:     Objective coefficient
%       If included, makes sure biomass is excluded from exch_rxns
%   model.lb:    Lower bounds (if revFlag=0)
%       Must be included if would like to know med_rxns
%   model.ub:    Lower bounds (if revFlag=1)
% revFlag: Indicates if the exchange reactions are reversed ("1") or not ("0")
%
%OUTPUT
% exch_rxns: Vector with indices of the exchange reactions
% med_rxns: Vector with indices of the medium reactions (model can use)
%   Need to include model.lb/model.ub to compute
%
% Meghan Thommes 07/14/2017
% Meghan Thommes 03/03/2015

%% Check Inputs

if (nargin < 1)
    error('myfuns:findExchRxns:NotEnoughInputs', ...
        'Not enough inputs: need a model file');
elseif (nargin == 1)
    if ~isstruct(model)
        error('myfuns:findExchRxns:IncorrectInput', ...
            '"model" needs to be a structure');
    elseif ~isfield(model,'S')
        error('myfuns:findExchRxns:IncorrectInput', ...
            '"model" needs "S"');
    end
    revFlag = 0;
else
    if ~isstruct(model)
        error('myfuns:findExchRxns:IncorrectInput', ...
            '"model" needs to be a structure');
    elseif ~isfield(model,'S')
        error('myfuns:findExchRxns:IncorrectInput', ...
            '"model" needs "S" field');
    end
end

if (nargin > 1)
    if ~isscalar(revFlag)
        error('myfuns:findExchRxns:IncorrectInput', ...
            '"revFlag" needs to be a scalar');
    elseif (revFlag ~= 0 && revFlag ~= 1)
        error('myfuns:findExchRxns:IncorrectInput', ...
            '"revFlag" needs to be "0" or "1"');
    end
end

%% Find Exchange Reactions

if ~revFlag % exchange reactions are not reversed
    % Find Reactions (Columns) with Stoichiometry of -1
    [~,C_1] = find(model.S == -1);
    C_1 = unique(C_1);
else % exchange reactions are reversed
    % Find Reactions (Columns) with Stoichiometry of 1
    [~,C_1] = find(model.S == 1);
    C_1 = unique(C_1);
end
% Count How Many Non-Zero Elements Reactions Have
[~,C_0] = find(model.S);
C = unique(C_0);
N = histcounts(C_0,C);
% Reactions with only ONE Non-Zero Elements that are Equal to +/- 1 is an Exchange Reaction
exch_rxns = intersect(C_1,C(N == 1));

% Eliminate Objective Function
if isfield(model, 'c')
    exch_rxns(intersect(find(model.c),exch_rxns)) = [];
end

% Find Medium Reactions
if ~revFlag % exchange reactions are not reversed
    % Find Reactions with Lower Bound less than Zero
    if isfield(model, 'lb')
        med_rxns = exch_rxns(model.lb(exch_rxns) < 0);
    else
        med_rxns = [];
    end
else % exchange reactions are reversed
    % Find Reactions with Upper Bound greater than Zero
    if isfield(model, 'ub')
        med_rxns = exch_rxns(model.ub(exch_rxns) > 0);
    else
        med_rxns = [];
    end
end

end