function model = sbml2struct(fileName)
%%sbml2struct Convert SBML into a MATLAB structure
% model = sbml2struct(fileName)
%
%REQUIRED INPUT
% filename            File name for SBML to convert
%
%OUTPUT
% Returns a MATLAB structure with the following fields:
% model.description   Model name (string)
% model.ID            Model ID (string)
% model.S             Stoichiometric matrix (m x n)
% model.mets          Metabolite IDs (m x 1)
% model.metNames      Metabolite names (m x 1)
% model.metCompts     Metabolite compartments (m x 1)
% model.metFormulas   Metabolite chemical formula (m x 1)
% model.rxns          Reaction IDs (n x 1)
% model.rxnNames      Reaction names (n x 1)
% model.rxnCompts     Reaction compartments (n x 1)
% model.lb            Lower bounds (n x 1)
% model.ub            Upper bounds (n x 1)
% model.rev           Reaction reversibility (n x 1)
% model.c             Objective coefficients: 1 if objective, 0 if other (n x 1)
%
% 01/04/2016 Meghan Thommes - added ID field
% 09/15/2015 Meghan Thommes

%% Validate Inputs

% Check that file exists
if ~exist(fileName,'file')
    error('myFuns:sbml2struct:IncorrectType', ...
        ['Input file ' fileName ' not found']);
end

%%

MATLAB_model = sbml2model(TranslateSBML(fileName));

model.description = MATLAB_model.name; % model name
model.ID = MATLAB_model.ID; % model ID
model.S = MATLAB_model.S; % stoichiometric matrix (m x n)

[numMets,numRxns] = size(model.S);

% Metabolites
model.mets = cell(numMets,1);
model.metNames = cell(numMets,1);
model.metCompts = cell(numMets,1);
for met = 1:numMets
    model.mets{met} = MATLAB_model.Metabolites(met).ID; % metabolite IDs (m x 1)
    model.metNames{met} = MATLAB_model.Metabolites(met).Name; % metabolite names (m x 1)
    model.metCompts{met} = MATLAB_model.Metabolites(met).Compartment; % metabolite compartments (m x 1)
    model.metFormulas{met} = MATLAB_model.Metabolites(met).Formula;   % metabolite chemical formula (m x 1)
end

% Reactions
model.rxns = cell(numRxns,1);
model.rxnNames = cell(numRxns,1);
model.rxnCompts = cell(numRxns,1);
model.lb = zeros(numRxns,1);
model.ub = zeros(numRxns,1);
model.rev = zeros(numRxns,1);
model.c = zeros(numRxns,1);
for rxn = 1:numRxns
    model.rxns{rxn} = MATLAB_model.Reactions(rxn).ID; % reaction IDs (n x 1)
    model.rxnNames{rxn} = MATLAB_model.Reactions(rxn).Name; % eaction names (n x 1)
    model.rxnCompts{rxn} = MATLAB_model.Reactions(rxn).Compartment; % reaction compartments (n x 1)
    model.lb(rxn) = MATLAB_model.Reactions(rxn).LowerBound; % lower bounds (n x 1)
    model.ub(rxn) = MATLAB_model.Reactions(rxn).UpperBound; % upper bounds (n x 1)
    model.rev(rxn) = MATLAB_model.Reactions(rxn).Reversible; % reaction reversibility (n x 1)
    model.c(rxn) = MATLAB_model.Reactions(rxn).ObjectiveCoefficient; % objective coefficients (n x 1)
end


end
