function [MURI_models,base_MURI_model] = make_MURI_models(model,data,rxn_names,rxn_bounds)
%MAKE_MURI_MODELS Create metabolic models based on the binary variables
%from the MURI data
%   [MURI_models,base_MURI_model] = make_MURI_models(model,data)
%   [MURI_models,base_MURI_model] = make_MURI_models(model,data,rxn_names,rxn_bounds)
%
%REQUIRED INPUTS
% model: metabolic model used in the MURI algorithm to find binary
%   variables. Must contain the following fields:
%       S:  stoichiometric matrix
%       c: objective coefficients
%       mets: metabolite names (short)
%       metNames: metabolite names (long)
%       rxns: reaction names (short)
%       rxnNames: reaction names (long)
%   	lb: lower bounds
%   	ub: upper bounds
% data: Structure from the output of process_flux_all. Must contain the
%   following fields:
%       sparse_con: sparsity constraints, vector
%       mets: metabolite names, vector
%       rxns: reaction names, vector
%       biomass: biomass flux, vector
%       flux: reaction fluxes, matrix [rxns x sparse_con]
%       int: reaction binary variables t, matrix [rxns x sparse_con]
%       exch_idx: index of exchange ("uptake") reactions
%       trspt_idx: index of transport ("exchange") reactions
%       intl_idx: index of intracellular ("internal") reactions
%
%OPTIONAL INPUTS
% rxn_names: name of reactions that need their lower and/or upper bounds
%   changed
% rxn_bounds: matrix of new lower and upper bounds [lb, ub]
%
%OUTPUT
% MURI_models: cell structure of the metabolic models for each of the
%   sparsity constraints in data.sparse_con. Contains the following fields:
%       S:  stoichiometric matrix
%       c: objective coefficients
%       mets: metabolite names (short)
%       metNames: metabolite names (long)
%       rxns: reaction names (short)
%       rxnNames: reaction names (long)
%   	lb: lower bounds
%   	ub: upper bounds
% base_MURI_model: metabolic model without any of the sparsity constraints.
%   Contains the following fields:
%       S:  stoichiometric matrix
%       c: objective coefficients
%       mets: metabolite names (short)
%       metNames: metabolite names (long)
%       rxns: reaction names (short)
%       rxnNames: reaction names (long)
%   	lb: lower bounds
%   	ub: upper bounds
%       exch_idx: index of exchange ("uptake") reactions
%       trspt_idx: index of transport ("exchange") reactions
%       intl_idx: index of intracellular ("internal") reactions
%
% Meghan Thommes 03/29/2017

%% Base Model

tol = 1E-6;

% Stoichiometric Matrix
base_MURI_model.S = model.S;
% Metabolites
[~,idx_core_mets,idx_data_mets] = intersect(model.mets,data.mets,'stable'); % keep in same order as E. coli core model
base_MURI_model.mets = data.mets(idx_data_mets);
base_MURI_model.metNames = model.metNames(idx_core_mets);
% Reactions
[~,idx_core_rxns,idx_data_rxns] = intersect(model.rxns,data.rxns,'stable'); % keep in same order as E. coli core model
base_MURI_model.rxns = data.rxns(idx_data_rxns);
base_MURI_model.rxnNames = model.rxnNames(idx_core_rxns);
% RHS (dx/dt)
base_MURI_model.b = zeros(size(base_MURI_model.mets));
% Upper & Lower Bounds
base_MURI_model.lb = model.lb(idx_core_rxns);
base_MURI_model.ub = model.ub(idx_core_rxns);
% Objective
base_MURI_model.c = model.c(idx_core_rxns);
% Check Model
sol_base_model = FBA(base_MURI_model,'',true);
sol_model = FBA(model,'',true);
if abs(sol_model.objectiveValue - sol_base_model.objectiveValue) > tol
    disp('Error in Base Model Creation: Objective value does not match model')
end
% Reaction Indices
exch_idx = zeros(size(base_MURI_model.rxns)); exch_idx(data.exch_idx) = 1; base_MURI_model.exch_idx = find(exch_idx(idx_data_rxns)==1);
trspt_idx = zeros(size(base_MURI_model.rxns)); trspt_idx(data.trspt_idx) = 1; base_MURI_model.trspt_idx = find(trspt_idx(idx_data_rxns)==1);
intl_idx = zeros(size(base_MURI_model.rxns)); intl_idx(data.intl_idx) = 1; base_MURI_model.intl_idx = find(intl_idx(idx_data_rxns)==1);

% Update Exchange Lower & Upper Bounds
if exist('rxn_names')
    for ii = 1:numel(rxn_names)
        [~,idx_rxn,~] = intersect(base_MURI_model.rxns,rxn_names{ii});
        base_MURI_model.lb(idx_rxn) = rxn_bounds(ii,1);
        base_MURI_model.ub(idx_rxn) = rxn_bounds(ii,2);
    end
end

%% Model with Sparsity Constraints

for ii = 1:numel(data.sparse_con)
    MURI_models{ii} = base_MURI_model;   
    MURI_models{ii}.lb = data.int(idx_data_rxns,ii).*base_MURI_model.lb;
    MURI_models{ii}.ub = data.int(idx_data_rxns,ii).*base_MURI_model.ub;
end

end

