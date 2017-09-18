function [model_flux] = processAlgData(alg_data,model,mediumRxns,tol)
%%processAlgData Reformat algorithm output
%
% model_flux = processAlgData(alg_data,model)
% model_flux = processAlgData(alg_data,model,mediumRxns,tol)
%
%REQUIRED INPUTS
% alg_data: Algorithm output. Structure must contain fields:
%   mets: metabolite names (short) [cell array]
%   rxns: reaction names (short) [cell array]
%   sparse_con: intracellular sparsity constraint [vector]
%   trspt_con: transport sparsity constraint [vector, size of sparse_con]
%   opt_status: secondary optimization flag (1 if min sum of fluxes, 0 if infeasible) [vector, size of sparse_con]
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
% Meghan Thommes 09/05/2017

%% Check Inputs and Assign Variables

if (nargin < 2)
    error('myfuns:processAlgData:NotEnoughInputs', ...
        'Not enough inputs: need a "alg_data" & "model"');
end
if ~isstruct(alg_data)
    error('myfuns:processAlgData:IncorrectType', ...
            '"alg_data" needs to be a structure');
elseif ~isfield(alg_data,'mets') || ~isfield(alg_data,'rxns') || ~isfield(alg_data,'sparse_con') || ~isfield(alg_data,'trspt_con') || ~isfield(alg_data,'opt_status') || ~isfield(alg_data,'model')
    error('myfuns:processAlgData:IncorrectType', ...
        '"alg_data" needs "mets", "rxns", "sparse_con", "trspt_con", "opt_status", & "model" fields');
end

numModels = numel(alg_data.model);

for model_num = 1:numModels
    if ~isstruct(alg_data.model{model_num})
        error('myfuns:processAlgData:IncorrectType', ...
            ['"alg_data.model{' int2str(model_num) '}" needs to be a structure']);
    elseif ~isfield(alg_data.model{model_num},'biomass') || ~isfield(alg_data.model{model_num},'flux') || ~isfield(alg_data.model{model_num},'int')
        error('myfuns:processAlgData:IncorrectType', ...
        ['"alg_data.model{' int2str(model_num) '}" needs "biomass", "flux", & "int" fields']);
    end
end

% Medium
if ~exist('mediumRxns','var')
    [~,exch_rxns] = identifyExchRxns(model);
    [exch_mets,] = identifyExchMets(model,exch_rxns);
    mediumMets = model.mets(exch_mets);
    mediumRxns = model.rxns(exch_rxns);
else
    [~,~,medRxns_idx] = intersect(mediumRxns,model.rxns,'stable');
    if numel(medRxns_idx) ~= numel(mediumRxns)
        error('myfuns:processAlgData:IncorrectInput',...
            'Error in mediumRxns: Not all reactions are in model')
    end
    [exch_mets,~] = identifyExchMets(model,medRxns_idx);
    mediumMets = model.mets(exch_mets);
end

% Tolerance from Zero (+/-)
if ~exist('tol','var')
    tol = 1E-6;
end

%%

% Indices of Exchange, Transport, and Intracellular Reactions
exch_idx = identifyExchRxns(model);
trspt_idx = identifyTrsptRxns(model);
intl_idx = setdiff(setdiff(1:numel(alg_data.rxns),exch_idx),trspt_idx); if(size(intl_idx,2) > size(intl_idx,1)); intl_idx = intl_idx'; end
% Check that Have Indices for All Reactions
if numel(alg_data.rxns) ~= (numel(exch_idx) + numel(trspt_idx) + numel(intl_idx))
    error('myfuns:algorithm2models:IncorrectCalc', ...
        'Incorrect calculation of reaction indices');
end
% Remove Biomass Reaction from Intracellular Reactions Indices
[~,rm_idx,~] = intersect(intl_idx,alg_data.biomass_id);
intl_idx(rm_idx) = [];

% Metabolite Names
[~,~,metsIdx_model] = intersect(alg_data.mets,model.mets,'stable'); % keep in the same order as the algorithm data
alg_data.metNames = model.metNames(metsIdx_model);
if numel(metsIdx_model) ~= numel(model.mets)
    error('myfuns:algorithm2models:IncorrectCalc', ...
        'Incorrect calculation of metabolite indices');
end

% Reaction Names
[~,~,rxnsIdx_model] = intersect(alg_data.rxns,model.rxns,'stable'); % keep in the same order as the algorithm data
alg_data.rxnNames = model.rxnNames(rxnsIdx_model);
if numel(rxnsIdx_model) ~= numel(model.rxns)
    error('myfuns:algorithm2models:IncorrectCalc', ...
        'Incorrect calculation of reaction indices');
end

%% Reorder Flux and Integer Vectors

% Pre-Allocate
model_flux = cell(numModels,1);

% Flux and Integer
totalCalcExchFlux = zeros(numel(exch_idx),numel(alg_data.sparse_con));
for model_num = 1:numModels
    if ~issorted(alg_data.sparse_con) % if not listed in ascending order, change to ascending order
        model_flux{model_num}.intl_con = fliplr(alg_data.sparse_con);
        model_flux{model_num}.trspt_con = fliplr(alg_data.trspt_con);
        model_flux{model_num}.opt_status = fliplr(alg_data.opt_status);
        model_flux{model_num}.mets = alg_data.mets;
        model_flux{model_num}.metNames = alg_data.metNames;
        model_flux{model_num}.rxns = alg_data.rxns;
        model_flux{model_num}.rxnNames = alg_data.rxnNames;
        model_flux{model_num}.biomass = fliplr(alg_data.biomass(model_num,:));
        model_flux{model_num}.flux = fliplr(alg_data.model{model_num}.flux);
        model_flux{model_num}.flux(abs(model_flux{model_num}.flux) < tol) = 0;
        model_flux{model_num}.int = fliplr(alg_data.model{model_num}.int);
        model_flux{model_num}.int = model_flux{model_num}.int > tol;
        model_flux{model_num}.exch_idx = exch_idx;
        model_flux{model_num}.trspt_idx = trspt_idx;
        model_flux{model_num}.intl_idx = intl_idx;
        model_flux{model_num}.mediumMets = mediumMets;
        model_flux{model_num}.mediumRxns = mediumRxns;
        model_flux{model_num}.warningFlag.BinaryVariableFlux = zeros(1,numel(alg_data.sparse_con));
        model_flux{model_num}.warningFlag.CalcExchFlux = zeros(1,numel(alg_data.sparse_con));
    else % if not listed in ascending order, change to ascending order
        model_flux{model_num}.intl_con = alg_data.sparse_con;
        model_flux{model_num}.trspt_con = alg_data.trspt_con;
        model_flux{model_num}.opt_status = alg_data.opt_status;
        model_flux{model_num}.mets = alg_data.mets;
        model_flux{model_num}.metNames = alg_data.metNames;
        model_flux{model_num}.rxns = alg_data.rxns;
        model_flux{model_num}.rxnNames = alg_data.rxnNames;
        model_flux{model_num}.biomass = alg_data.biomass(model_num,:);
        model_flux{model_num}.flux = alg_data.model{model_num}.flux;
        model_flux{model_num}.flux(abs(model_flux{model_num}.flux) < tol) = 0;
        model_flux{model_num}.int = alg_data.model{model_num}.int;
        model_flux{model_num}.int = model_flux{model_num}.int > tol;
        model_flux{model_num}.exch_idx = exch_idx;
        model_flux{model_num}.trspt_idx = trspt_idx;
        model_flux{model_num}.intl_idx = intl_idx;
        model_flux{model_num}.mediumMets = mediumMets;
        model_flux{model_num}.mediumRxns = mediumRxns;
        model_flux{model_num}.warningFlag.BinaryVariableFlux = zeros(1,numel(alg_data.sparse_con));
        model_flux{model_num}.warningFlag.CalcExchFlux = zeros(1,numel(alg_data.sparse_con));
    end
    
    % Set Flux and Binary Variables to Zero if Biomass is Below a Certain Threshold
    if ~isempty(intersect(find(model_flux{model_num}.biomass < tol  | model_flux{model_num}.biomass < tol),find(model_flux{model_num}.biomass ~= 0)))
        intlCon_bioIdx = intersect(find(model_flux{model_num}.biomass < tol  | model_flux{model_num}.biomass < tol),find(model_flux{model_num}.biomass ~= 0));
        model_flux{model_num}.biomass(intlCon_bioIdx) = 0;
        model_flux{model_num}.flux(:,intlCon_bioIdx) = 0;
        model_flux{model_num}.int(:,intlCon_bioIdx) = 0;
    end
    
    % Add No Growth Case
    if model_flux{model_num}.biomass(1)~=0
        model_flux{model_num}.intl_con = [model_flux{model_num}.intl_con(1)-1, model_flux{model_num}.intl_con];
        model_flux{model_num}.trspt_con = [model_flux{model_num}.trspt_con(1), model_flux{model_num}.trspt_con];
        model_flux{model_num}.opt_status = [model_flux{model_num}.opt_status(1), model_flux{model_num}.opt_status];
        model_flux{model_num}.biomass = [0, model_flux{model_num}.biomass];
        model_flux{model_num}.flux = [zeros(numel(model_flux{model_num}.rxns),1), model_flux{model_num}.flux];
        model_flux{model_num}.int = [zeros(numel(model_flux{model_num}.rxns),1), model_flux{model_num}.int];
        model_flux{model_num}.warningFlag.BinaryVariableFlux = [0, model_flux{model_num}.warningFlag.BinaryVariableFlux];
        model_flux{model_num}.warningFlag.CalcExchFlux = [0, model_flux{model_num}.warningFlag.CalcExchFlux];
        if model_num==1
            totalCalcExchFlux = [zeros(numel(exch_idx),1), totalCalcExchFlux];
        end
    end
    
    % Check that Flux=0 when t_i=0
    t0_flux_idx = find(abs(model_flux{model_num}.flux(model_flux{model_num}.int==0)) > tol);
    if ~isempty(t0_flux_idx)
        warning('myfuns:algorithm2models:ConstraintViolation', ...
            ['Model ' int2str(model_num) ' has non-zero flux when t_i=0 -- Updated flux to zero'])
        temp_flux = zeros(size(model_flux{model_num}.flux));
        temp_flux(model_flux{model_num}.int==0) = model_flux{model_num}.flux(model_flux{model_num}.int==0);
        model_flux{model_num}.flux(find(temp_flux)) = 0;
        [rxns_idx,intlCon_idx] = find(temp_flux); intlCon_idx = unique(intlCon_idx); rxns_idx = unique(rxns_idx);
        model_flux{model_num}.warningFlag.BinaryVariableFlux(intlCon_idx) = 1;
    end
    
    % Calculate Exchange Flux from Transport Flux
    [exchMets_idx,~] = identifyExchMets(model,model_flux{model_num}.exch_idx);
    trsptS = model.S(exchMets_idx,model_flux{model_num}.trspt_idx);
    calcExchFlux = trsptS*model_flux{model_num}.flux(model_flux{model_num}.trspt_idx,:); % exchange flux = stoichiometric matrix * transport flux
    calcExchFlux(abs(calcExchFlux) < tol) = 0;
    totalCalcExchFlux = totalCalcExchFlux + calcExchFlux;
    
    % Update Exchange Flux
    model_flux{model_num}.totalExchFlux = model_flux{model_num}.flux(model_flux{model_num}.exch_idx,:);
    model_flux{model_num}.flux(model_flux{model_num}.exch_idx,:) = calcExchFlux;
end

% Check Calculated Exchange Fluxes
diff_flux = model_flux{1}.totalExchFlux - totalCalcExchFlux;
tempS = zeros(size(trsptS)); tempS(find(trsptS)) = 1;
if ~isempty(find(abs(diff_flux) > sum(tempS,2).*tol)) % error can accumlate with each additional transport reactions involved, acceptable error is within #trspt_rxns*tolerance
    warning('myfuns:algorithm2models:IncorrectCalc', ...
        ['Incorrect calculation of exchange flux. Max difference = ' num2str(max(abs(diff_flux(:))))]);
    [~,intlCon_idx] = find(abs(diff_flux) > sum(tempS,2).*tol);
    for model_num = 1:numModels
        model_flux{model_num}.warningFlag.CalcExchFlux(intlCon_idx) = 1;
    end
end

end

