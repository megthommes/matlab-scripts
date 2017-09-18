function [exch_rxns,med_rxns] = identifyExchRxns(model,exch_mets)
%%identifyExchRxns Identify exchange and "medium" reactions.
%   [exch_rxns,med_rxns] = identifyExchRxns(model)
%   [exch_rxns,med_rxns] = identifyExchRxns(model,exch_mets)
%
%By default, identifyExchRxns assumes that exchange reactions are written
%as export reactions: 'A <==>'; therefore, uptaking a metabolite is a
%negative flux and the stoichiometric coefficient of metabolite A is -1.
%
%REQUIRED INPUT
% The model structure must contain the following fields:
%   model.S: Stoichiometric matrix
%
%OPTIONAL INPUTS
% The model structure my contain the following fields:
%   model.lb: Lower bounds. Must be included if would like to know med_rxns
%   model.c: Objective coefficient. If included, makes sure biomass is 
%       excluded from exch_rxns
% exch_mets: Exchange metabolite indices of interest (numeric array)
%
%OUTPUT
% exch_rxns: Vector with indices of the exchange reactions
% med_rxns: Vector with indices of the medium reactions (model can use)
%   Need to include model.lb to compute
%
% Meghan Thommes 07/14/2017

%% Check Inputs

if (nargin < 1)
    error('myfuns:identifyExchRxns:NotEnoughInputs', ...
        'Not enough inputs: need a model file');
elseif (nargin >= 1)
    if ~isstruct(model)
        error('myfuns:identifyExchRxns:IncorrectInput', ...
            '"model" needs to be a structure');
    elseif ~isfield(model,'S')
        error('myfuns:identifyExchRxns:IncorrectInput', ...
            '"model" needs "S"');
    end
end

%% Identify Exchange Reactions

% Find Reactions (Columns) with Stoichiometry of -1
[~,C_1] = find(model.S == -1);
C_1 = unique(C_1);
% Count How Many Non-Zero Elements Reactions Have
[~,C_0] = find(model.S);
C = unique(C_0);
N = histcounts(C_0,C);
% Reactions with only ONE Non-Zero Elements that are Equal to - 1 are Exchange Reactions
exch_rxns = intersect(C_1,C(N == 1));

% Eliminate Objective Function
if isfield(model, 'c')
    exch_rxns(intersect(find(model.c),exch_rxns)) = [];
end

% Find Medium Reactions
if isfield(model, 'lb')
    % Find Reactions with Lower Bound less than Zero
    med_rxns = exch_rxns(model.lb(exch_rxns) < 0);
else
    med_rxns = [];
end

if exist('exch_mets','var')
    clear C_1
    % Find Reactions (Columns) with Stoichiometry of -1
    [~,C_1] = find(model.S(exch_mets,:) == -1);
    C_1 = unique(C_1);
    exch_rxns = intersect(C_1,exch_rxns);
    
    if numel(exch_rxns) ~= numel(exch_mets)
        error('myfuns:identifyExchRxns:IncorrectCalc', ...
            'Number of exch_rxns not equal to number of exch_mets');
    end
end

end