function [model_flux] = identifyExchangedMets(model_flux,tol)
%%identifyExchangedMets Identify exchange metabolites
%
% model_flux = identifyExchangedMets(model_flux,tol)
% model_flux = identifyExchangedMets(model_flux)
%
%REQUIRED INPUTS
% model_flux: Cell structure, size of the number of models reactions are allocated
%  to. Contains fields:
%   intl_con: intracellular sparsity constraint [vector]
%   trspt_con: transport sparsity constraint [vector, size of sparse_con]
%   opt_status: secondary optimization flag (1 if min sum of fluxes, 0 if infeasible) [vector, size of sparse_con]
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
%OPTIONAL INPUT
% tol: tolerance from zero (default = 1E-6)
%
%OUTPUT
% model_flux: Added exchMets structure with fields:
%   metNames: metabolite names [cell array] (metabolites x 1)
%   fluxType: Number indicates what type of flux occurs: [matrix] (metabolites x intracellular sparsity constraint)
%       0: not used by model 1 or model 2 (->)
%       1: produced by model 1 and model 2 (1,2->)
%       2: produced by model 1 and not used by model 2 (1->)
%       3: not used by model 1 and produced by model 2 (2->)
%       4: consumed by model 1 and model 2 (->1,2)
%       5: consumed by model 1 and not used by model 2 (->1)
%       6: not used by model 1 and consumed by model 2 (->2)
%       7: consumed by model 1 and produced by model 2 (2->1)
%       8: produced by model 1 and consumed by model 2 (1->2)
%   model1_to_model2: 1 if exchanged, 0 if not exchanged [matrix] (metabolites x intracellular sparsity constraint)
%   model2_to_model1: 1 if exchanged, 0 if not exchanged [matrix] (metabolites x intracellular sparsity constraint)
%
% Meghan Thommes 09/07/2017

%% Check Inputs and Assign Variables

if (nargin < 2)
    error('myfuns:processAlgData:NotEnoughInputs', ...
        'Not enough inputs: need a "model_flux" & "mediumRxns"');
end
for model_num = 1:2
    if ~isstruct(model_flux{model_num})
        error('myfuns:reassignModels:IncorrectType', ...
            ['"model_flux{' int2str(model_num) '}" needs to be a structure']);
    elseif ~isfield(model_flux{model_num},'biomass') || ~isfield(model_flux{model_num},'flux') || ~isfield(model_flux{model_num},'int')
        error('myfuns:reassignModels:IncorrectType', ...
        ['"model_flux{' int2str(model_num) '}" needs "biomass", "flux", & "int" fields']);
    end
end

% Tolerance from Zero (+/-)
if ~exist('tol','var')
    tol = 1E-6;
end

% Medium
mediumRxns = model_flux{1}.mediumRxns;

%% Classify Exchanged Flux

fluxType = classifyExchangeFlux(model_flux{1}.flux(model_flux{1}.exch_idx,:),model_flux{2}.flux(model_flux{2}.exch_idx,:),tol);

%% Find Exchanged Metabolites

% Metabolites Secreted by Model 2 and Used by Model 1
[idx_metNames21,idx_intlCon21] = find(fluxType == 7);
rm_idx = find(ismember(model_flux{1}.rxns(model_flux{1}.exch_idx(idx_metNames21)), mediumRxns));
idx_metNames21(rm_idx) = []; idx_intlCon21(rm_idx) = []; % remove metabolites that are in the medium
exchMets_21 = zeros(size(fluxType));
for ii = 1:numel(idx_metNames21)
    exchMets_21(idx_metNames21(ii),idx_intlCon21(ii)) = 1;
end

% Metabolites Secreted by Model 1 and Used by Model 2
[idx_metNames12,idx_intlCon12] = find(fluxType == 8);
rm_idx = find(ismember(model_flux{1}.rxns(model_flux{1}.exch_idx(idx_metNames12)), mediumRxns));
idx_metNames12(rm_idx) = []; idx_intlCon12(rm_idx) = []; % remove metabolites that are in the medium
exchMets_12 = zeros(size(fluxType));
for ii = 1:numel(idx_metNames12)
    exchMets_12(idx_metNames12(ii),idx_intlCon12(ii)) = 1;
end

%% Update model_flux

for model_num = 1:2
    model_flux{model_num}.exchMets.metNames = model_flux{model_num}.metNames(model_flux{model_num}.exch_idx);
    model_flux{model_num}.exchMets.fluxType = fluxType;
    model_flux{model_num}.exchMets.model1_to_model2 = exchMets_12;
    model_flux{model_num}.exchMets.model2_to_model1 = exchMets_21;
end

end