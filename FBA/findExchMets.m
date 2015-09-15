function [exch_mets,nut_mets] = findExchMets(model,revFlag)
%findExchMets Find exchange and "nutrient" metabolites
%[exch_mets,nut_mets] = findExchMets(model,revFlag)
%
%By default, findExchMets assumes that exchange reactions are written as
%export reactions: 'A <==>'; therefore, uptaking a metabolite is a
%negative flux and the stoichiometric coefficient of metabolite A is -1.
%
%REQUIRED INPUT
% The model structure must contain the following fields:
%   model.S:     Stoichiometric matrix
%   model.c:     Objective coefficients
%   model.lb:    Lower bounds
%
%OPTIONAL INPUT
% revFlag indicates if the exchange reactions are reversed ("1") or not ("0")
%
%OUTPUT
% exch_mets is a Boolean vector indicating whether each metabolite in the
% model is an exchange ("1") or not ("0")
%   default: "0" (false)
%
% nut_mets is a Boolean vector indicating whether each metabolite in the
% model is a nutrient ("1") or not ("0")
%   Nutrient: Lower bound is less than zero
%
%3/4/2015 - Meghan Thommes

%%

% Check Inputs
if (nargin < 1)
    error('Not enough inputs: need a model file');
elseif (nargin == 1)
    if ~isstruct(model)
        error('"model" needs to be a structure');
    elseif ~isfield(model,'S') || ~isfield(model,'c') || ~isfield(model,'lb')
        error('"model" needs "S", "c", and "lb" fields');
    end
    revFlag = 0;
else
    if ~isstruct(model)
        error('"model" needs to be a structure');
    elseif ~isfield(model,'S') || ~isfield(model,'c') || ~isfield(model,'lb')
        error('"model" needs "S", "c", and "lb" fields');
    end
end

if (nargin > 1)
    if ~isscalar(revFlag)
        error('"revFlag" needs to be a scalar');
    elseif (revFlag ~= 0 && revFlag ~= 1)
        error('"revFlag" needs to be "0" or "1"');
    end
end

% Find exchange reactions
[exch_rxns,nut_rxns] = findExchRxns(model,revFlag);

mets_matrix = model.S(:,exch_rxns);
[exch_mets,~] = find(mets_matrix ~= 0);

nuts_matrix = model.S(:,nut_rxns);
[nut_mets,~] = find(nuts_matrix ~= 0);

end

