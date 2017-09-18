function [exch_mets,med_mets] = identifyExchMets(model,exch_rxns)
%identifyExchMets Identify exchange and "medium" metabolites.
%   [exch_mets,med_mets] = identifyExchMets(model)
%   [exch_mets,med_mets] = identifyExchMets(model,exch_rxns)
%
%By default, identifyExchMets assumes that exchange reactions are written as
%export reactions: 'A <==>'; therefore, uptaking a metabolite is a
%negative flux and the stoichiometric coefficient of metabolite A is -1.
%
%REQUIRED INPUT
% The model structure must contain the following field:
%   model.S:     Stoichiometric matrix
%
%OPTIONAL INPUTS
% The model structure my contain the following fields:
%   model.lb: Lower bounds. Must be included if would like to know med_rxns
%   model.c: Objective coefficient. If included, makes sure biomass is 
%       excluded from exch_rxns
% exch_rxns: Exchange reaction indices of interest (numeric array)
%
%OUTPUT
% exch_mets: Vector with indices of the exchange metabolites
% med_rxns: Vector with indices of the medium metabolites (model can use)
%   Need to include model.lb to compute
%
% Meghan Thommes 08/14/2017

%% Check Inputs

if (nargin < 1)
    error('myfuns:identifyExchMets:NotEnoughInputs', ...
        'Not enough inputs: need a model file');
elseif (nargin >= 1)
    if ~isstruct(model)
        error('myfuns:identifyExchMets:IncorrectInput', ...
            '"model" needs to be a structure');
    elseif ~isfield(model,'S')
        error('myfuns:identifyExchMets:IncorrectInput', ...
            '"model" needs "S"');
    end
end

%% Identify Exchange Metabolites

if ~exist('exch_rxns','var')
    % Identify Exchange Reactions
    [exch_rxns,med_rxns] = identifyExchRxns(model);
    
    % Identify Metabolites (Rows) with Non-Zero Stoichiometry
    mets_matrix = model.S(:,exch_rxns);
    [exch_mets,~] = find(mets_matrix);
    
    % Identify Medium Metabolites
    med_matrix = model.S(:,med_rxns);
    [med_mets,~] = find(med_matrix);
else
    % Identify Exchange Metabolites
    mets_matrix = model.S(:,exch_rxns);
    [exch_mets,~] = find(mets_matrix);
    
    if numel(exch_mets) ~= numel(exch_rxns)
        error('myfuns:identifyExchMets:IncorrectCalc', ...
            'Number of exch_mets not equal to number of exch_rxns');
    end
    
    % Identify Medium Metabolites
    [~,med_rxns] = identifyExchRxns(model);
    med_matrix = model.S(:,med_rxns);
    [med_mets,~] = find(med_matrix);
end

end

