function [dolmn_models,dolmn_base_model,model_flux] = algorithm2models(alg_data,model,model_description,mediumRxns,atpm_name)
%%algorithm2models Convert algorithm output to metabolic models
%
% [dolmn_models,dolmn_base_model,model_flux] = algorithm2models(alg_data,model,model_description)
% [dolmn_models,dolmn_base_model,model_flux] = algorithm2models(alg_data,model,model_description,mediumRxns,atpm_name)
%
%REQUIRED INPUTS
% alg_data: Algorithm output. Structure must contain fields:
%   mets: metabolite names (short) [cell array]
%   rxns: reaction names (short) [cell array]
%   sparse_con: intracellular sparsity constraint [vector]
%   trspt_con: transport sparsity constraint [vector, size of sparse_con]
%   model: cell structure, size of the number of models reactions are allocated
%    to. Must contain fields:
%       biomass: biomass flux [vector, size of sparse_con]
%       flux: reaction flux [matrix, size of rxns x size of sparse_con]
%       int: reaction binary variables [matrix, size of rxns x size of sparse_con]
% model: Metabolic model used in the algorithm for reaction allocation.
%  Must contain the fields:
%   S: Stoichiometric matrix [sparse matrix, number of metabolites x number of reactions]
%   mets: metabolite names (short) [cell array]
%   metNames: metabolite names (long) [cell array]
%   rxns: reaction names (short) [cell array]
%   rxnNames: reaction names (long) [cell array]
%   c: objective coefficients [vector, length of rxns]
%   lb: reaction lower bounds [vector, length of rxns]
%   ub: reaction upper bounds [vector, length of rxns]
% model_description: description of model (e.g. E. coli core 1 model)
%
%OPTIONAL INPUTS
% mediumRxns: Names of reactions in the medium (default are the
%  exchange reactions with a non-zero lower bound in model)
%   - Must correspond to model.rxns
% atpm_name: name of the ATPM reaction (default = 'ATPM')
%
%OUTPUTS
% dolmn_models: cell structure of the metabolic models for each of the
% sparsity constraints in model_flux.intl_con. # of models x 1 cell array.
% Contains the following fields:
%   warningFlag: 0 if no warning, 1 if warning
%      FBA_sol: The model does not grow
%      Sv_viol: Constraint violation, S*v != 0
%      bounds_viol: Constraint violation, lb !<= v !<= ub
%   S:  stoichiometric matrix
%   mets: metabolite names (short)
%   metNames: metabolite names (long)
%   rxns: reaction names (short)
%   rxnNames: reaction names (long)
%   b: right hand side
%   lb: lower bounds
%   ub: upper bounds
%   c: objective coefficients
%   exch_idx: index of exchange reactions
%   trspt_idx: index of transport reactions
%   intl_idx: index of intracellular reactions
%   intl_con: intracellular sparsity constraint
%   trspt_con: transport sparsity constraint
%   flux: reaction fluxes from dolmn
%   int: reaction binary variables from dolmn
%   sol: objective value from FBA
%   exchMets_AB/exchMets_BA: exchanged metabolites between model A and model B
%       **Only if there are two or more models!
% dolmn_base_model: metabolic model without any of the sparsity constraints.
%   Contains the following fields:
%   S:  stoichiometric matrix
%   c: objective coefficients
%   mets: metabolite names (short)
%   metNames: metabolite names (long)
%   rxns: reaction names (short)
%   rxnNames: reaction names (long)
%   lb: lower bounds
%   ub: upper bounds
%   exch_idx: index of exchange reactions
%   trspt_idx: index of transport reactions
%   intl_idx: index of intracellular reactions
%   warningFlag: 0 if no warning, 1 if warning
%      FBA_sol: The objective value does not match model
% model_flux: Cell structure, size of the number of models reactions are allocated
%  to. Contains fields:
%   intl_con: intracellular sparsity constraint [vector]
%   trspt_con: transport sparsity constraint [vector, size of sparse_con]
%   mets: metabolite names (short) [cell array]
%   metNames: metabolite names (long) [cell array]
%   rxns: reaction names (short) [cell array]
%   rxnNames: reaction names (long) [cell array]
%   biomass: biomass flux [vector, size of sparse_con]
%   flux: reaction flux [matrix, size of rxns x size of sparse_con]
%   int: reaction binary variables [matrix, size of rxns x size of sparse_con]
%   mediumMets: Names of metabolites in the medium
%   mediumRxns: Names of reactions in the medium
%   warningFlag: 0 if no warning, 1 if warning. Structure with fields:
%      BinaryVariableFlux: Constraint violation, model has non-zero flux
%       when t_i=0 [vector] (size of intl_con)
%      CalcExchFlux: Incorrect calculation of exchange flux [cell array] (size of intl_con)
%      FBA_sol: The constrained model does not grow at that intracellular 
%       sparsity constraint [vector] (size of intl_con)
%      sol: The FBA solution of the constrained model does at that
%       intracellular sparsity constraint [vector] (size of intl_con)
%      Sv_viol: Constraint violation, S*v != 0, at that intracellular 
%       sparsity constraint [matrix] (number of metabolites x length of intl_con)
%      Sv: S*v at that intracellular sparsity constraint [matrix] (number of metabolites x length of intl_con)
%      bounds_viol: Constraint violation, lb !<= v !<= ub, at that 
%       intracellular sparsity constraint [vector] (size of intl_con)
%      intl_bounds: Which internal reaction bound is violated [matrix] (number of reactions x length of intl_con)
%   totalExchFlux: Original exchange flux from the community stoichiometric
%    matrix
%   exchMets: Cell structure, contains the following fields:
%       metNames: metabolite names [cell array] (metabolites x 1)
%       fluxType: Number indicates what type of flux occurs [matrix] (metabolites x intracellular sparsity constraint)
%           0: not used by model 1 or model 2 (->)
%           1: produced by model 1 and model 2 (1,2->)
%           2: produced by model 1 and not used by model 2 (1->)
%           3: not used by model 1 and produced by model 2 (2->)
%           4: consumed by model 1 and model 2 (->1,2)
%           5: consumed by model 1 and not used by model 2 (->1)
%           6: not used by model 1 and consumed by model 2 (->2)
%           7: consumed by model 1 and produced by model 2 (2->1)
%           8: produced by model 1 and consumed by model 2 (1->2)
%       modelA_to_modelB/modelB_to_modelA: metabolites exchanged between
%           model A and model B at each intracellular sparsity constraint.
%           1 if exchanged, 0 if not exchanged [matrix]
%           (metabolites x intracellular sparsity constraint)
%
% Meghan Thommes 09/05/2017

%% Variables

% Number of Models
numModels = numel(alg_data.model);

% Tolerance from Zero (+/-)
tol = 1E-4;

% Medium
if ~exist('mediumRxns','var')
    [~,exch_rxns] = identifyExchRxns(model);
    mediumRxns = model.rxns(exch_rxns);
else
    [~,~,medRxns_idx] = intersect(mediumRxns,model.rxns,'stable');
    if numel(medRxns_idx) ~= numel(mediumRxns)
        error('myfuns:algorithm2models:Incorrectinput',...
            'Error in mediumRxns: Not all reactions are in model')
    end
end

%% Add ATPM Reaction if Not Present

% ATPM Reaction
if ~exist('atpm_name','var')
    atpm_name = 'ATPM';
end
if isempty(intersect(alg_data.rxns,atpm_name))
    [~,atpm_idx,~] = intersect(model.rxns,atpm_name);
    atpm_value = model.lb(atpm_idx);
    if atpm_idx == 1
        for model_num = 1:numModels
            alg_data.rxns = [atpm_name; alg_data.model{model_num}.rxns];
            alg_data.model{model_num}.rxns = alg_data.rxns;
            alg_data.model{model_num}.flux = [atpm_value.*ones(1,numel(alg_data.sparse_con)); alg_data.model{model_num}.flux];
            alg_data.model{model_num}.flux(atpm_idx,alg_data.model{model_num}.biomass==0) = 0; % value = 0 if no biomass
            alg_data.model{model_num}.int = [ones(1,numel(alg_data.sparse_con)); alg_data.model{model_num}.int];
            alg_data.model{model_num}.int(atpm_idx,alg_data.model{model_num}.biomass==0) = 0; % value = 0 if no biomass
        end
    elseif atpm_idx > 1 && atpm_idx < numel(alg_data.rxns)
        for model_num = 1:numModels
            alg_data.rxns = [alg_data.model{model_num}.rxns(1:atpm_idx-1); atpm_name; ...
                alg_data.model{model_num}.rxns(atpm_idx:end)];
            alg_data.model{model_num}.rxns = alg_data.rxns;
            alg_data.model{model_num}.flux = [alg_data.model{model_num}.flux(1:atpm_idx-1,:); ...
                atpm_value.*ones(1,numel(alg_data.sparse_con)); alg_data.model{model_num}.flux(atpm_idx:end,:)];
            alg_data.model{model_num}.flux(atpm_idx,alg_data.model{model_num}.biomass==0) = 0; % value = 0 if no biomass
            alg_data.model{model_num}.int = [alg_data.model{model_num}.int(1:atpm_idx-1,:); ...
                ones(1,numel(alg_data.sparse_con)); alg_data.model{model_num}.int(atpm_idx:end,:)];
            alg_data.model{model_num}.int(atpm_idx,alg_data.model{model_num}.biomass==0) = 0; % value = 0 if no biomass
        end
    elseif atpm_idx == numel(alg_data.rxns)
        for model_num = 1:numModels
            alg_data.rxns = [alg_data.model{model_num}.rxns; atpm_name];
            alg_data.model{model_num}.rxns = alg_data.rxns;
            alg_data.model{model_num}.flux = [alg_data.model{model_num}.flux; atpm_value.*ones(1,numel(alg_data.sparse_con))];
            alg_data.model{model_num}.flux(atpm_idx,alg_data.model{model_num}.biomass==0) = 0; % value = 0 if no biomass
            alg_data.model{model_num}.int = [alg_data.model{model_num}.int; ones(1,numel(alg_data.sparse_con))];
            alg_data.model{model_num}.int(atpm_idx,alg_data.model{model_num}.biomass==0) = 0; % value = 0 if no biomass
        end
    end
    % Add 1 to Account for ATPM Reaction
    alg_data.sparse_con = alg_data.sparse_con + 1;
end

%% Process Algorithm Data

model_flux = processAlgData(alg_data,model,mediumRxns,tol);

%% Reassign Model Labels and Identify Exchanged Metabolites

if numModels == 2
    % Reassign Model Labels
    model_flux = reassignModels(model_flux,'Flux');
    
    % Identify Exchanged Metabolites
    model_flux = identifyExchangedMets(model_flux,tol);
    for model_num = 1:numModels
        exchMets = model_flux{model_num}.exchMets;
        exchMets.fluxType12 = exchMets.fluxType;
        exchMets = rmfield(exchMets,'fluxType');
        model_flux{model_num}.exchMets = orderfields(exchMets,{'metNames', 'fluxType12','model1_to_model2','model2_to_model1'});
    end
    
elseif numModels == 3
    tempModel_flux = cell(2,1);
    
    % Identify Exchanged Metabolites - Models 1 & 2
    tempModel_flux{1} = model_flux{1};
    tempModel_flux{2} = model_flux{2};
    tempModel_flux = identifyExchangedMets(tempModel_flux,tol);
    for model_num = 1:numModels
        exchMets = tempModel_flux{1}.exchMets;
        exchMets.fluxType12 = exchMets.fluxType;
        exchMets = rmfield(exchMets,'fluxType');
        exchMets.model1_to_model2_save = exchMets.model1_to_model2;
        exchMets.model2_to_model1_save = exchMets.model2_to_model1;
        model_flux{model_num}.exchMets = exchMets;
    end
    
    % Identify Exchanged Metabolites - Models 1 & 3
    tempModel_flux{1} = model_flux{1};
    tempModel_flux{2} = model_flux{3};
    tempModel_flux = identifyExchangedMets(tempModel_flux,tol);
    for model_num = 1:numModels
        exchMets = tempModel_flux{1}.exchMets;
        exchMets.fluxType13 = exchMets.fluxType;
        exchMets = rmfield(exchMets,'fluxType');
        exchMets.model1_to_model3 = exchMets.model1_to_model2;
        exchMets = rmfield(exchMets,'model1_to_model2');
        exchMets.model3_to_model1 = exchMets.model2_to_model1;
        exchMets = rmfield(exchMets,'model2_to_model1');
        model_flux{model_num}.exchMets = exchMets;
    end
    
    % Identify Exchanged Metabolites - Models 2 & 3
    tempModel_flux{1} = model_flux{2};
    tempModel_flux{2} = model_flux{3};
    tempModel_flux = identifyExchangedMets(tempModel_flux,tol);
    for model_num = 1:numModels
        exchMets = tempModel_flux{1}.exchMets;
        exchMets.fluxType23 = exchMets.fluxType;
        exchMets = rmfield(exchMets,'fluxType');
        exchMets.model2_to_model3 = exchMets.model1_to_model2;
        exchMets = rmfield(exchMets,'model1_to_model2');
        exchMets.model3_to_model2 = exchMets.model2_to_model1;
        exchMets = rmfield(exchMets,'model2_to_model1');
        model_flux{model_num}.exchMets = exchMets;
    end
    
    for model_num = 1:numModels
        exchMets = model_flux{model_num}.exchMets;
        exchMets.model1_to_model2 = exchMets.model1_to_model2_save;
        exchMets = rmfield(exchMets,'model1_to_model2_save');
        exchMets.model2_to_model1 = exchMets.model2_to_model1_save;
        exchMets = rmfield(exchMets,'model2_to_model1_save');
        model_flux{model_num}.exchMets = orderfields(exchMets,{'metNames', 'fluxType12','model1_to_model2','model2_to_model1', 'fluxType13','model1_to_model3','model3_to_model1', 'fluxType23','model2_to_model3','model3_to_model2'});
    end
    
end

%% Make FBA Models

[dolmn_models,dolmn_base_model,model_flux] = makeFBAmodels(model,model_flux,model_description,mediumRxns,tol);
   
end