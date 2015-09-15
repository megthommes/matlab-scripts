function [exch_rxns,nut_rxns] = findExchRxns(model,revFlag)
%%findExchRxns Find exchange and "nutrient" reactions
%[exch_rxns,nut_rxns] = findExchRxns(model)
%
%By default, findExchRxns assumes that exchange reactions are written as
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
% exch_rxns is a Boolean vector indicating whether each reaction in the
% model is an exchange ("1") or not ("0")
%   default: "0" (false)
%
% nut_rxns is a Boolean vector indicating whether each reaction in the
% model is a nutrient ("1") or not ("0")
%   Nutrient: Lower bound is less than zero
%
%3/3/2015 - Meghan Thommes

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

% Examine stoichiometric matrix
if ~revFlag % exchange reactions are not reversed
    exch_rxns = full((sum(model.S==-1,1) == 1) & (sum(model.S~=0) == 1))';
    % (sum(model.S==-1,1) == 1): find columns (reactions) with stoichiometry of -1
    % (sum(model.S~=0,1) == 1): make sure columns (reactions) have only one non-zero element (metabolite)
else % exchange reactions are reversed
    exch_rxns = full((sum(model.S==1,1) == 1) & (sum(model.S~=0) == 1))';
    % (sum(model.S==1,1) == 1): find columns (reactions) with stoichiometry of 1
    % (sum(model.S~=0,1) == 1): make sure columns (reactions) have only one non-zero element (metabolite)
end

% Eliminate objective function
exch_rxns(model.c~=0) = false;

% Find "nutrient" reactions
nut_rxns = full(model.lb < 0 & exch_rxns);

end