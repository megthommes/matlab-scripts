function [model_flux] = reassignModels(model_flux,distance_metric)
%%reassignModels Reassign model flux and binary variables based on a
%%distance metric
%
% model_flux = reassignModels(model_flux)
% model_flux = reassignModels(model_flux,distance_metric)
%
%REQUIRED INPUT
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
%   mediumMets: Names of metabolites in the medium
%   mediumRxns: Names of reactions in the medium
%   warningFlag: 0 if no warning, 1 if warning. Structure with fields:
%      BinaryVariableFlux: Constraint violation, model has non-zero flux when t_i=0
%      CalcExchFlux: Incorrect calculation of exchange flux
%   totalExchFlux: Original exchange flux from the community stoichiometric
%    matrix
%
%OPTIONAL INPUT
% distance_metric: Metric to classify models by. [String] Can be:
%   Flux: Flux
%   Int: Binary variables (default)
%   IntFlux: Binarized flux (0 if zero flux, 1 if non-zero flux)
%
%OUTPUT
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
%   mediumRxns: Names of metabolites in the medium
%   mediumRxns: Names of reactions in the medium
%   warningFlag: 0 if no warning, 1 if warning. Structure with fields:
%      BinaryVariableFlux: Constraint violation, model has non-zero flux when t_i=0
%      CalcExchFlux: Incorrect calculation of exchange flux
%   totalExchFlux: Original exchange flux from the community stoichiometric
%    matrix
%
% Meghan Thommes 09/05/2017

%% Check Inputs and Assign Variables

if (nargin < 1)
    error('myfuns:reassignModels:NotEnoughInputs', ...
        'Not enough inputs: need a "model_flux"');
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

if ~exist('distance_metric','var')
    distance_metric = 'int';
else
    distance_metric = lower(distance_metric);
end

% Temporary Variables
M1_biomass = model_flux{1}.biomass;
M2_biomass = model_flux{2}.biomass;
M1_flux = model_flux{1}.flux;
M2_flux = model_flux{2}.flux;
M1_int = model_flux{1}.int;
M2_int = model_flux{2}.int;

%% Classify Based on Minimum Distance

if strcmp(distance_metric,'flux')
    for intlCon_num = find(M1_biomass,1):numel(M1_biomass)-1 % start at first index with biomass flux
        % Calculate Distance
        D11 = pdist2(M1_flux(:,intlCon_num)', M1_flux(:,intlCon_num+1)', 'euclidean'); % model 1 at this intl_con to the next intl_con
        D22 = pdist2(M2_flux(:,intlCon_num)', M2_flux(:,intlCon_num+1)', 'euclidean'); % model 2 at this intl_con to the next intl_con
        D12 = pdist2(M1_flux(:,intlCon_num)', M2_flux(:,intlCon_num+1)', 'euclidean'); % model 1 at this intl_con to model 2 at the next intl_con
        D21 = pdist2(M2_flux(:,intlCon_num)', M1_flux(:,intlCon_num+1)', 'euclidean'); % model 2 at this intl_con to model 1 at the next intl_con
        
        % Model Prediction
        if D12 < D11 || D21 < D22 % swap
            M1_biomass(:,intlCon_num+1) = model_flux{2}.biomass(intlCon_num+1);
            M2_biomass(:,intlCon_num+1) = model_flux{1}.biomass(intlCon_num+1);
            M1_flux(:,intlCon_num+1) = model_flux{2}.flux(:,intlCon_num+1);
            M2_flux(:,intlCon_num+1) = model_flux{1}.flux(:,intlCon_num+1);
            M1_int(:,intlCon_num+1) = model_flux{2}.int(:,intlCon_num+1);
            M2_int(:,intlCon_num+1) = model_flux{1}.int(:,intlCon_num+1);
        else % keep the same
            M1_biomass(:,intlCon_num+1) = model_flux{1}.biomass(intlCon_num+1);
            M2_biomass(:,intlCon_num+1) = model_flux{2}.biomass(intlCon_num+1);
            M1_flux(:,intlCon_num+1) = model_flux{1}.flux(:,intlCon_num+1);
            M2_flux(:,intlCon_num+1) = model_flux{2}.flux(:,intlCon_num+1);
            M1_int(:,intlCon_num+1) = model_flux{1}.int(:,intlCon_num+1);
            M2_int(:,intlCon_num+1) = model_flux{2}.int(:,intlCon_num+1);
        end
    end
elseif strcmp(distance_metric,'int')
    M1_int = double(M1_int);
    M2_int = double(M2_int);
    
    for intlCon_num = find(M1_biomass,1):numel(M1_biomass)-1 % start at first index with biomass flux
        % Calculate Distance
        D11 = pdist2(M1_int(:,intlCon_num)', M1_int(:,intlCon_num+1)', 'hamming'); % model 1 at this intl_con to the next intl_con
        D22 = pdist2(M2_int(:,intlCon_num)', M2_int(:,intlCon_num+1)', 'hamming'); % model 2 at this intl_con to the next intl_con
        D12 = pdist2(M1_int(:,intlCon_num)', M2_int(:,intlCon_num+1)', 'hamming'); % model 1 at this intl_con to model 2 at the next intl_con
        D21 = pdist2(M2_int(:,intlCon_num)', M1_int(:,intlCon_num+1)', 'hamming'); % model 2 at this intl_con to model 1 at the next intl_con
        
        % Model Prediction
        if D12 < D11 || D21 < D22 % swap
            M1_biomass(:,intlCon_num+1) = model_flux{2}.biomass(intlCon_num+1);
            M2_biomass(:,intlCon_num+1) = model_flux{1}.biomass(intlCon_num+1);
            M1_flux(:,intlCon_num+1) = model_flux{2}.flux(:,intlCon_num+1);
            M2_flux(:,intlCon_num+1) = model_flux{1}.flux(:,intlCon_num+1);
            M1_int(:,intlCon_num+1) = double(model_flux{2}.int(:,intlCon_num+1));
            M2_int(:,intlCon_num+1) = double(model_flux{1}.int(:,intlCon_num+1));
        else % keep the same
            M1_biomass(:,intlCon_num+1) = model_flux{1}.biomass(intlCon_num+1);
            M2_biomass(:,intlCon_num+1) = model_flux{2}.biomass(intlCon_num+1);
            M1_flux(:,intlCon_num+1) = model_flux{1}.flux(:,intlCon_num+1);
            M2_flux(:,intlCon_num+1) = model_flux{2}.flux(:,intlCon_num+1);
            M1_int(:,intlCon_num+1) = double(model_flux{1}.int(:,intlCon_num+1));
            M2_int(:,intlCon_num+1) = double(model_flux{2}.int(:,intlCon_num+1));
        end
    end
    
    M1_int = logical(M1_int);
    M2_int = logical(M2_int);
elseif strcmp(distance_metric,'intflux')
    % Binarize Flux
    M1_intFlux = zeros(size(M1_flux));
    M1_intFlux(M1_flux ~= 0) = 1;
    M2_intFlux = zeros(size(M2_flux));
    M2_intFlux(M2_flux ~= 0) = 1;
    
    for intlCon_num = find(M1_biomass,1):numel(M1_biomass)-1 % start at first index with biomass flux
        % Calculate Distance
        D11 = pdist2(M1_intFlux(:,intlCon_num)', M1_intFlux(:,intlCon_num+1)', 'hamming'); % model 1 at this intl_con to the next intl_con
        D22 = pdist2(M2_intFlux(:,intlCon_num)', M2_intFlux(:,intlCon_num+1)', 'hamming'); % model 2 at this intl_con to the next intl_con
        D12 = pdist2(M1_intFlux(:,intlCon_num)', M2_intFlux(:,intlCon_num+1)', 'hamming'); % model 1 at this intl_con to model 2 at the next intl_con
        D21 = pdist2(M2_intFlux(:,intlCon_num)', M1_intFlux(:,intlCon_num+1)', 'hamming'); % model 2 at this intl_con to model 1 at the next intl_con
        
        % Model Prediction
        if D12 < D11 || D21 < D22 % swap
            M1_biomass(:,intlCon_num+1) = model_flux{2}.biomass(intlCon_num+1);
            M2_biomass(:,intlCon_num+1) = model_flux{1}.biomass(intlCon_num+1);
            M1_flux(:,intlCon_num+1) = model_flux{2}.flux(:,intlCon_num+1);
            M2_flux(:,intlCon_num+1) = model_flux{1}.flux(:,intlCon_num+1);
            M1_int(:,intlCon_num+1) = model_flux{2}.int(:,intlCon_num+1);
            M2_int(:,intlCon_num+1) = model_flux{1}.int(:,intlCon_num+1);
            M1_intFlux(:,intlCon_num+1) = zeros(size(M1_flux(:,intlCon_num+1))); M1_intFlux(M1_flux(:,intlCon_num+1) ~= 0) = 1;
            M2_intFlux(:,intlCon_num+1) = zeros(size(M2_flux(:,intlCon_num+1))); M2_intFlux(M2_flux(:,intlCon_num+1) ~= 0) = 1;
        else % keep the same
            M1_biomass(:,intlCon_num+1) = model_flux{1}.biomass(intlCon_num+1);
            M2_biomass(:,intlCon_num+1) = model_flux{2}.biomass(intlCon_num+1);
            M1_flux(:,intlCon_num+1) = model_flux{1}.flux(:,intlCon_num+1);
            M2_flux(:,intlCon_num+1) = model_flux{2}.flux(:,intlCon_num+1);
            M1_int(:,intlCon_num+1) = model_flux{1}.int(:,intlCon_num+1);
            M2_int(:,intlCon_num+1) = model_flux{2}.int(:,intlCon_num+1);
            M1_intFlux(:,intlCon_num+1) = zeros(size(M1_flux(:,intlCon_num+1))); M1_intFlux(M1_flux(:,intlCon_num+1) ~= 0) = 1;
            M2_intFlux(:,intlCon_num+1) = zeros(size(M2_flux(:,intlCon_num+1))); M2_intFlux(M2_flux(:,intlCon_num+1) ~= 0) = 1;
        end
    end
end

% Model 1
model_flux{1}.biomass = M1_biomass;
model_flux{1}.flux = M1_flux;
model_flux{1}.int = M1_int;

% Model 2
model_flux{2}.biomass = M2_biomass;
model_flux{2}.flux = M2_flux;
model_flux{2}.int = M2_int;

end