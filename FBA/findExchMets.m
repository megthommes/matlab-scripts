function [exch_mets,med_mets] = findExchMets(model,revFlag)
%findExchMets Find exchange and metabolites and identify "medium" metabolites
%   [exch_mets,med_mets] = findExchMets(model)
%   [exch_mets,med_mets] = findExchMets(model,revFlag)
%
%By default, findExchMets assumes that exchange reactions are written as
%export reactions: 'A <==>'; therefore, uptaking a metabolite is a
%negative flux and the stoichiometric coefficient of metabolite A is -1.
%   If this is not the case, set revFlag=1.
%
%REQUIRED INPUT
% The model structure must contain the following field:
%   model.S:     Stoichiometric matrix
%
%OPTIONAL INPUTS
% The model structure my contain the following fields:
%   model.c:     Objective coefficient
%       If included, makes sure biomass is excluded from exch_rxns
%   model.lb:    Lower bounds (if revFlag=0)
%       Must be included if would like to know med_rxns
%   model.ub:    Lower bounds (if revFlag=1)
% revFlag: Indicates if the exchange reactions are reversed ("1") or not ("0")
%
%OUTPUT
% exch_mets: Vector with indices of the exchange metabolites
% med_rxns: Vector with indices of the medium metabolites (model can use)
%   Need to include model.lb/model.ub to compute
%
% Meghan Thommes 07/14/2017
% Meghan Thommes 03/04/2015

%% Check Inputs

if (nargin < 1)
    error('myfuns:findExchMets:NotEnoughInputs', ...
        'Not enough inputs: need a model file');
elseif (nargin == 1)
    if ~isstruct(model)
        error('myfuns:findExchMets:IncorrectInput', ...
            '"model" needs to be a structure');
    elseif ~isfield(model,'S')
        error('myfuns:findExchMets:IncorrectInput', ...
            '"model" needs "S"');
    end
    revFlag = 0;
else
    if ~isstruct(model)
        error('myfuns:findExchMets:IncorrectInput', ...
            '"model" needs to be a structure');
    elseif ~isfield(model,'S')
        error('myfuns:findExchMets:IncorrectInput', ...
            '"model" needs "S" field');
    end
end

if (nargin > 1)
    if ~isscalar(revFlag)
        error('myfuns:findExchMets:IncorrectInput', ...
            '"revFlag" needs to be a scalar');
    elseif (revFlag ~= 0 && revFlag ~= 1)
        error('myfuns:findExchMets:IncorrectInput', ...
            '"revFlag" needs to be "0" or "1"');
    end
end

%% Find Exchange Metabolites

% Find Exchange Reactions
[exch_rxns,med_rxns] = findExchRxns(model,revFlag);

% Find Metabolites (Rows) with Non-Zero Stoichiometry
mets_matrix = model.S(:,exch_rxns);
[exch_mets,~] = find(mets_matrix);

% Find Medium Metabolites
med_matrix = model.S(:,med_rxns);
[med_mets,~] = find(med_matrix);

end

