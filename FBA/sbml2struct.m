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
% model.b             Right-hand side of FBA problem (m x 1)
% model.rxns          Reaction IDs (n x 1)
% model.rxnNames      Reaction names (n x 1)
% model.rxnCompts     Reaction compartments (n x 1)
% model.lb            Lower bounds (n x 1)
% model.ub            Upper bounds (n x 1)
% model.rev           Reaction reversibility (n x 1)
% model.c             Objective coefficients: 1 if objective, 0 if other (n x 1)
%
% 01/08/2016 Meghan Thommes - eliminated loops
% 01/04/2016 Meghan Thommes - added ID & b fields
% 09/15/2015 Meghan Thommes

%% Validate Inputs

% Check that file exists
if ~exist(fileName,'file')
    error('myFuns:sbml2struct:IncorrectType', ...
        ['Input file ' fileName ' not found']);
end

%%

MATLAB_model = sbml2model(TranslateSBML(fileName));

model.description = MATLAB_model.Name; % model name
model.ID = MATLAB_model.ID; % model ID
model.S = MATLAB_model.S; % stoichiometric matrix (m x n)

[numMets,numRxns] = size(model.S);

% Metabolites
model.mets = {MATLAB_model.Metabolites.ID}'; % metabolite IDs (m x 1)
model.metNames = {MATLAB_model.Metabolites.Name}'; % metabolite names (m x 1)
model.metCompts = {MATLAB_model.Metabolites.Compartment}'; % metabolite compartments (m x 1)
model.metFormulas = {MATLAB_model.Metabolites.Formula}'; % metabolite chemical formula (m x 1)
model.b = zeros(numMets,1);

% Reactions
model.rxns = {MATLAB_model.Reactions.ID}'; % reaction IDs (n x 1)
model.rxnNames = {MATLAB_model.Reactions.Name}'; % reaction names (n x 1)
model.rxnCompts = {MATLAB_model.Reactions.Compartment}'; % reaction compartments (n x 1)
model.lb = [MATLAB_model.Reactions.LowerBound]'; % lower bounds (n x 1)
model.ub = [MATLAB_model.Reactions.UpperBound]'; % upper bounds (n x 1)
model.rev = [MATLAB_model.Reactions.Reversible]'; % reaction reversibility (n x 1)
model.c = [MATLAB_model.Reactions.ObjectiveCoefficient]'; % objective coefficients (n x 1)

end
