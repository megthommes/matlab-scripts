function [cmda_models,cmda_base_model,model_flux] = makeFBAmodels(model,model_flux,model_description,mediumRxns,tol)
%%makeFBAmodels Make FBA models from algorithm output.
%
% [cmda_models,cmda_base_model,model_flux] = makeFBAmodels(model,model_flux,model_description)
% [cmda_models,cmda_base_model,model_flux] = makeFBAmodels(model,model_flux,model_description,mediumRxns,tol)
%
%REQUIRED INPUTS
% alg_data: Algorithm output. Structure must contain fields:
%   mets: metabolite names (short) [cell array]
%   rxns: reaction names (short) [cell array]
%   sparse_con: intracellular sparsity constraint [vector]
%   trspt_con: transport sparsity constraint [vector, size of sparse_con]
%   bio_lb: lower bound on biomass flux
%   model: Cell structure, size of the number of models reactions are allocated
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
%
%OPTIONAL INPUTS
% mediumRxns: Names of reactions in the medium (default are the
%  exchange reactions with a non-zero lower bound in model)
%   - Must correspond to model.rxns
% tol: tolerance from zero (default = 1E-6)
%
%OUTPUTS
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
%   exch_idx: exchange reaction indices (vector)
%   trspt_idx: transport reaction indices (vector)
%   intl_idx: intracellular reaction indices (vector)
%   mediumMets: Names of metabolites in the medium
%   mediumRxns: Names of reactions in the medium
%   warningFlag: 0 if no warning, 1 if warning. Structure with fields:
%      BinaryVariableFlux: Constraint violation, model has non-zero flux when t_i=0
%      CalcExchFlux: Incorrect calculation of exchange flux
%   totalExchFlux: Original exchange flux from the community stoichiometric
%    matrix
%
% Meghan Thommes

%% Check Inputs and Assign Variables

numModels = numel(model_flux);

% Pre-Allocate
cmda_base_model = cell(numModels,1);
cmda_models = cell(numModels,1);

% Warning Flag: 1 if the objective value does not match model
for model_num = 1:numModels
    cmda_base_model{model_num}.warningFlag.FBA_sol = 0;
end

% Feasibility Tolerance
FBA_params.FeasibilityTol = 1e-2; % (max value)

%% Base Model

for model_num = 1:numModels
    % Metabolites
    cmda_base_model{model_num}.mets = model_flux{model_num}.mets;
    cmda_base_model{model_num}.metNames = model_flux{model_num}.metNames;
    % Reactions
    cmda_base_model{model_num}.rxns = model_flux{model_num}.rxns;
    cmda_base_model{model_num}.rxnNames = model_flux{model_num}.rxnNames;
    % Pathways
    if isfield(model,'subSystems')
        [~,~,rxnsIdx_model] = intersect(model_flux{model_num}.rxns,model.rxns,'stable'); % keep in the same order as the algorithm data
        cmda_base_model{model_num}.subSystems = model.subSystems(rxnsIdx_model);
    end
    % Stoichiometric Matrix
    cmda_base_model{model_num}.S = model.S;
    % Upper & Lower Bounds
    cmda_base_model{model_num}.lb = model.lb;
    cmda_base_model{model_num}.ub = model.ub;
    % RHS (dx/dt)
    cmda_base_model{model_num}.b = zeros(size(cmda_base_model{model_num}.mets));
    % Objective
    cmda_base_model{model_num}.c = model.c;
    % Reversibility Flag
    if isfield(model,'rev')
        [~,~,rxnsIdx_model] = intersect(model_flux{model_num}.rxns,model.rxns,'stable'); % keep in the same order as the algorithm data
        cmda_base_model{model_num}.rev = model.rev(rxnsIdx_model);
    end
    % Model Descripton
    cmda_base_model{model_num}.description = model_description;
    % Check Model
    sol_baseModel = FBA(cmda_base_model{model_num},'',true);
    sol_model = FBA(model,'',true);
    if abs(sol_model.objectiveValue - sol_baseModel.objectiveValue) > tol
        warning('myfuns:make_FBA_models:IncorrectCalc',...
            ['Error in Base Model ' int2str(model_num) ' Creation: Objective value does not match model'])
        cmda_base_model{model_num}.warningFlag.FBA_sol = 1;
    end
    
    %% Model with Sparsity Constraints
    
    cmda_base_model{model_num}.exch_idx = identifyExchRxns(cmda_base_model{model_num});
    cmda_base_model{model_num}.trspt_idx = identifyTrsptRxns(cmda_base_model{model_num});
    cmda_base_model{model_num}.intl_idx = setdiff(setdiff(1:numel(cmda_base_model{model_num}.rxns),cmda_base_model{model_num}.exch_idx),cmda_base_model{model_num}.trspt_idx); if size(cmda_base_model{model_num}.intl_idx,2) > size(cmda_base_model{model_num}.intl_idx,1); cmda_base_model{model_num}.intl_idx = cmda_base_model{model_num}.intl_idx'; end
    cmda_base_model{model_num}.bio_idx = find(cmda_base_model{model_num}.c);
    % Check that Have Indices for All Reactions
    if numel(cmda_base_model{model_num}.rxns) ~= (numel(cmda_base_model{model_num}.exch_idx) + numel(cmda_base_model{model_num}.trspt_idx) + numel(cmda_base_model{model_num}.intl_idx))
        error('myfuns:algorithm2models:IncorrectCalc', ...
            'Incorrect calculation of reaction indices');
    end
    % Remove Biomass Reaction from Intracellular Reactions Indices
    [~,rm_idx,~] = intersect(cmda_base_model{model_num}.intl_idx,cmda_base_model{model_num}.bio_idx);
    cmda_base_model{model_num}.intl_idx(rm_idx) = []; % remove biomass
    
    % Pre-Allocate
    intlCon_idx = find(model_flux{model_num}.biomass); % only make models when algorithm finds growth
    model_flux{model_num}.warningFlag.FBA_sol = zeros(1,numel(model_flux{model_num}.intl_con));
    model_flux{model_num}.warningFlag.sol = zeros(1,numel(model_flux{model_num}.intl_con));
    model_flux{model_num}.warningFlag.Sv_viol = zeros(1,numel(model_flux{model_num}.intl_con));
    model_flux{model_num}.warningFlag.Sv = zeros(numel(model_flux{model_num}.mets),numel(model_flux{model_num}.intl_con));
    model_flux{model_num}.warningFlag.bounds_viol = zeros(1,numel(model_flux{model_num}.intl_con));
    model_flux{model_num}.warningFlag.intl_bounds = zeros(numel(model_flux{model_num}.intl_idx),numel(model_flux{model_num}.intl_con));
    model_flux{model_num}.warningFlag.rev_viol = zeros(1,numel(model_flux{model_num}.intl_con));
    cmda_models{model_num} = cell(numel(intlCon_idx),1);
    for intlCon_num = 1:numel(intlCon_idx)
        cmda_models{model_num}{intlCon_num} = cmda_base_model{model_num};
        cmda_models{model_num}{intlCon_num}.intl_con = model_flux{model_num}.intl_con(intlCon_idx(intlCon_num));
        cmda_models{model_num}{intlCon_num}.trspt_con = model_flux{model_num}.trspt_con(intlCon_idx(intlCon_num));
        cmda_models{model_num}{intlCon_num}.flux = model_flux{model_num}.flux(:,intlCon_idx(intlCon_num));
        cmda_models{model_num}{intlCon_num}.int = model_flux{model_num}.int(:,intlCon_idx(intlCon_num));
        cmda_models{model_num}{intlCon_num}.exch_idx = cmda_base_model{model_num}.exch_idx;
        cmda_models{model_num}{intlCon_num}.trspt_idx = cmda_base_model{model_num}.trspt_idx;
        cmda_models{model_num}{intlCon_num}.intl_idx = cmda_base_model{model_num}.intl_idx;
        % If t_i=0, set lb=ub=0. If t_i=1, do not change lb or ub.
        cmda_models{model_num}{intlCon_num}.lb = cmda_models{model_num}{intlCon_num}.int.*cmda_base_model{model_num}.lb;
        cmda_models{model_num}{intlCon_num}.ub = cmda_models{model_num}{intlCon_num}.int.*cmda_base_model{model_num}.ub;
        % Add Information on Exchanged Metabolites to FBA Models
        if numModels == 2
            % Metabolite Names
            [exchMets_idx,~] = identifyExchMets(model,model_flux{model_num}.exch_idx);
            exchMetsNames = model_flux{model_num}.metNames(exchMets_idx);
            % Model 1 to Model 2
            [idx_exchMets,~] = find(model_flux{model_num}.exchMets.model1_to_model2); idx_exchMets = unique(idx_exchMets);
            cmda_models{model_num}{intlCon_num}.exchMets_12 = exchMetsNames(idx_exchMets);
            % Model 2 to Model 1
            [idx_exchMets,~] = find(model_flux{model_num}.exchMets.model2_to_model1); idx_exchMets = unique(idx_exchMets);
            cmda_models{model_num}{intlCon_num}.exchMets_21 = exchMetsNames(idx_exchMets);
        elseif numModels == 3
            % Metabolite Names
            [exchMets_idx,~] = identifyExchMets(model,model_flux{model_num}.exch_idx);
            exchMetsNames = model_flux{model_num}.metNames(exchMets_idx);
            % Model 1 to Model 2
            [idx_exchMets,~] = find(model_flux{model_num}.exchMets.model1_to_model2); idx_exchMets = unique(idx_exchMets);
            cmda_models{model_num}{intlCon_num}.exchMets_12 = exchMetsNames(idx_exchMets);
            % Model 2 to Model 1
            [idx_exchMets,~] = find(model_flux{model_num}.exchMets.model2_to_model1); idx_exchMets = unique(idx_exchMets);
            cmda_models{model_num}{intlCon_num}.exchMets_21 = exchMetsNames(idx_exchMets);
            % Model 1 to Model 3
            [idx_exchMets,~] = find(model_flux{model_num}.exchMets.model1_to_model3); idx_exchMets = unique(idx_exchMets);
            cmda_models{model_num}{intlCon_num}.exchMets_13 = exchMetsNames(idx_exchMets);
            % Model 3 to Model 1
            [idx_exchMets,~] = find(model_flux{model_num}.exchMets.model3_to_model1); idx_exchMets = unique(idx_exchMets);
            cmda_models{model_num}{intlCon_num}.exchMets_31 = exchMetsNames(idx_exchMets);
            % Model 2 to Model 3
            [idx_exchMets,~] = find(model_flux{model_num}.exchMets.model2_to_model3); idx_exchMets = unique(idx_exchMets);
            cmda_models{model_num}{intlCon_num}.exchMets_23 = exchMetsNames(idx_exchMets);
            % Model 3 to Model 2
            [idx_exchMets,~] = find(model_flux{model_num}.exchMets.model3_to_model2); idx_exchMets = unique(idx_exchMets);
            cmda_models{model_num}{intlCon_num}.exchMets_32 = exchMetsNames(idx_exchMets);
        end
        
        % Check if Constraints are Violated
        cmda_models{model_num}{intlCon_num}.warningFlag.FBA_sol = 0;
        cmda_models{model_num}{intlCon_num}.warningFlag.Sv_viol = 0;
        cmda_models{model_num}{intlCon_num}.warningFlag.bounds_viol = 0;
        cmda_models{model_num}{intlCon_num}.warningFlag.rev_viol = 0;
        % S*v=0
        Sv = cmda_models{model_num}{intlCon_num}.S * cmda_models{model_num}{intlCon_num}.flux;
        model_flux{model_num}.warningFlag.Sv(:,intlCon_idx(intlCon_num)) = Sv;
        Sv_viol = find(abs(Sv) > 5*tol);
        if ~isempty(Sv_viol)
            warning('myfuns:make_FBA_models:ConstraintViolation',...
                ['Model ' int2str(model_num) ' violates |S*v|<=' num2str(5*tol) ' for intlCon=' int2str(cmda_models{model_num}{intlCon_num}.intl_con) ' & trsptCon=' int2str(cmda_models{model_num}{intlCon_num}.trspt_con) '. Max Sv = ' num2str(max(abs(Sv(:))))]);
            cmda_models{model_num}{intlCon_num}.warningFlag.Sv_viol = 1;
            model_flux{model_num}.warningFlag.Sv_viol(intlCon_idx(intlCon_num)) = 1;
        end
        % lb <= v <= ub
        bounds_viol = find(cmda_models{model_num}{intlCon_num}.flux(cmda_models{model_num}{intlCon_num}.intl_idx) < cmda_models{model_num}{intlCon_num}.lb(cmda_models{model_num}{intlCon_num}.intl_idx) | ...
            cmda_models{model_num}{intlCon_num}.flux(cmda_models{model_num}{intlCon_num}.intl_idx) > cmda_models{model_num}{intlCon_num}.ub(cmda_models{model_num}{intlCon_num}.intl_idx));
        if ~isempty(bounds_viol)
            warning('myfuns:make_FBA_models:ConstraintViolation',...
                ['Model ' int2str(model_num) ' violates lb <= v <= ub for intlCon=' int2str(cmda_models{model_num}{intlCon_num}.intl_con) ' & trsptCon=' int2str(cmda_models{model_num}{intlCon_num}.trspt_con)]);
            cmda_models{model_num}{intlCon_num}.warningFlag.bounds_viol = 1;
            model_flux{model_num}.warningFlag.bounds_viol(intlCon_idx(intlCon_num)) = 1;
            model_flux{model_num}.warningFlag.intl_bounds(bounds_viol,intlCon_idx(intlCon_num)) = 1;
        end
        
        % Check if Model Grows
        exch_rxns = cmda_models{model_num}{intlCon_num}.exch_idx;
        [~,medium_idx,~] = intersect(cmda_models{model_num}{intlCon_num}.rxns,mediumRxns); % medium indices
        exchMets_idx = setdiff(exch_rxns(find(cmda_models{model_num}{intlCon_num}.flux(exch_rxns) < 0)),medium_idx); % metabolites taken up by the model but not in the medium
        [~,exchRxns_idx,~] = intersect(exch_rxns,exchMets_idx); % indices to exchange reactions
        exchRxns_lb = cmda_models{model_num}{intlCon_num}.lb(exch_rxns); % medium
        exchRxns_lb(exchRxns_idx) = cmda_models{model_num}{intlCon_num}.flux(exchMets_idx); % exchanged metabolites
        exchRxns_ub = cmda_models{model_num}{intlCon_num}.ub(exch_rxns);
        sol = FBA(cmda_models{model_num}{intlCon_num}, [exch_rxns exchRxns_lb exchRxns_ub], ...
            1, FBA_params);
        cmda_models{model_num}{intlCon_num}.sol = sol.objectiveValue;
        model_flux{model_num}.warningFlag.sol(intlCon_idx(intlCon_num)) = sol.objectiveValue;
        if sol.objectiveValue == 0 || isnan(sol.objectiveValue)
            % Decrease Lower Bound for Metabolites are Taken Up...
            exchRxns_lb(exchRxns_lb < 0) = 1.05.*exchRxns_lb(exchRxns_lb < 0);
            % ...Does the Model Grow Now?
            sol = FBA(cmda_models{model_num}{intlCon_num}, [exch_rxns exchRxns_lb exchRxns_ub], ...
                1, FBA_params);
            cmda_models{model_num}{intlCon_num}.sol = sol.objectiveValue;
            model_flux{model_num}.warningFlag.sol(intlCon_idx(intlCon_num)) = sol.objectiveValue;
            % If the Model Doesn't Grow Now, Flag It
            if sol.objectiveValue == 0 || isnan(sol.objectiveValue)
                warning('myfuns:make_FBA_models:ConstraintViolation',...
                    ['Model ' int2str(model_num) ' does not grow for intlCon=' int2str(cmda_models{model_num}{intlCon_num}.intl_con) ' & trsptCon=' int2str(cmda_models{model_num}{intlCon_num}.trspt_con)]);
                cmda_models{model_num}{intlCon_num}.warningFlag.FBA_sol = 1;
                model_flux{model_num}.warningFlag.FBA_sol(intlCon_idx(intlCon_num)) = 1;
                cmda_models{model_num}{intlCon_num}.warningFlag.sol = 0;
                model_flux{model_num}.warningFlag.sol(intlCon_idx(intlCon_num)) = 0;
            end
        end
    end
end