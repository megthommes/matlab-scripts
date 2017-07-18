function [model_flux] = process_flux_all_new(flux_all_data,atpm_info)
%%PROCESS_FLUX_ALL process flux_all into data that can be parsed
%
% model = PROCESS_FLUX_ALL(flux_all,sparse_con,S_info,rxns,numModels)
%
%REQUIRED INPUT
% flux_all_data: structure of Taiyao's output, contains the following fields:
%   mets: metabolite names, vector
%   rxns: reaction names, vector
%   sparse_con: sparsity constraints, vector
%   flux_all: output of the reaction-constraint algorithm, matrix [flux_all_row_id x sparse_con]
%   biomass: biomass flux, vector
%   biomass_id: index of the biomass reaction
%   flux_all_row_id: reaction type, # indicates model number, vector
%       U: exchange/uptake
%       E#: transport/exchange
%       I#: intracellular/internal
%       T#: reaction binary variables
%   model: cell structure that contains the following fields:
%       sparse_con: sparsity constraints, vector
%       rxns: reaction names, vector
%       biomass: biomass flux, vector
%       flux: reaction fluxes, matrix [rxns x sparse_con]
%       int: reaction binary variables t, matrix [rxns x sparse_con]
%   S: stoichiometric matrix
% atpm_info: structure, contains the following fields:
%       atpm_name: name of ATPM reaction
%       atpm_value: value of ATPM reaction
%
%OUTPUT
% model_flux: cell structure that contains the following fields:
%   sparse_con: sparsity constraints, vector
%   mets: metabolite names, vector
%   rxns: reaction names, vector
%   biomass: biomass flux, vector
%   flux: reaction fluxes, matrix [rxns x sparse_con]
%   int: reaction binary variables t, matrix [rxns x sparse_con]
%   exch_idx: index of exchange ("uptake") reactions
%   trspt_idx: index of transport ("exchange") reactions
%   intl_idx: index of intracellular ("internal") reactions
%
% Meghan Thommes 05/12/2017 - Calculate exchange flux from transport flux
%                             in 2 models
% Meghan Thommes 03/29/2017 - Updated int so that it is actually of the 
%                             size [rxns x sparse_con]
% Meghan Thommes 03/03/2017 - Exclude biomass reaction from intl_idx
% Meghan Thommes 02/28/2017

%% Reassign Model Labels

tol = 1E-3;
numModels = numel(flux_all_data.model);
temp_model = struct;
temp_model.S = flux_all_data.S;

% Pre-Allocate
biomass = cell(numModels,1);
flux = cell(numModels,1);
int = cell(numModels,1);

% Sparsity Constraints
exch_idx = findExchRxns(temp_model);
trspt_idx = findTrsptRxns(temp_model);
intl_idx = setdiff(setdiff(1:numel(flux_all_data.rxns),exch_idx),trspt_idx);
sparse_con = flux_all_data.sparse_con;
mets = flux_all_data.mets;
rxns = flux_all_data.rxns;
for ii = 1:numModels
    biomass{ii} = flux_all_data.biomass(ii,:);
    flux{ii} = flux_all_data.model{ii}.flux;
    flux{ii}(abs(flux{ii}) < tol) = 0;
    if size(flux_all_data.model{ii}.int,1) == size(flux{ii})
        int{ii} = flux_all_data.model{ii}.int;
    else
        int{ii} = [ones(numel(exch_idx)+numel(trspt_idx),numel(sparse_con)); flux_all_data.model{ii}.int];
    end
end
% Remove Biomass Reaction from Intracellular/Internal Reactions
[~,biomass_idx,~] = intersect(intl_idx,flux_all_data.biomass_id);
intl_idx(biomass_idx) = [];

% If not listed in ascending order, change to ascending order
if ~issorted(sparse_con) 
    sparse_con = fliplr(sparse_con);
    for ii = 1:numModels
        biomass{ii} = fliplr(biomass{ii});
        flux{ii} = fliplr(flux{ii});
        int{ii} = fliplr(int{ii});
    end
end

if numModels == 1 % 1 Model
    model_flux{1}.sparse_con = sparse_con;
    model_flux{1}.mets = mets;
    model_flux{1}.rxns = rxns;
    model_flux{1}.biomass = biomass{1};
    model_flux{1}.flux = flux{1};   
    model_flux{1}.int = int{1};
    model_flux{1}.exch_idx = exch_idx;
    model_flux{1}.trspt_idx = trspt_idx;
    model_flux{1}.intl_idx = intl_idx;
    
    [exch_flux1,exch_mets] = trsptFlux2exchFlux(mets,flux_all_data.S,trspt_idx,flux{1}(trspt_idx,:));
    exch_flux1(abs(exch_flux1) < tol) = 0;
    [~,exch_mets_idx,~] = intersect(mets,exch_mets);
    exchRxns_metsIdx = findExchRxnsFromExchMets(temp_model,exch_mets_idx);
    % Check new exchange fluxes
    total_flux = flux{1}(exchRxns_metsIdx,:);
    calc_flux = exch_flux1;
    diff_flux = total_flux - calc_flux;
    if ~isempty(find(abs(diff_flux) > tol))
        warning('myfuns:process_flux_all_new:IncorrectCalc', ...
            'Incorrect calculation of exchange flux');
    end
    model_flux{1}.total_flux = total_flux;
    model_flux{1}.calc_flux = calc_flux;
    
    % Add ATPM Reaction
    model_flux{1}.rxns{end+1} = atpm_info.atpm_name;
    model_flux{1}.flux(end+1,:) = atpm_info.atpm_value;
    model_flux{1}.flux(end,model_flux{1}.biomass==0) = 0; % value = 0 if no biomass
    model_flux{1}.int(end+1,:) = 1;
    model_flux{1}.int(end,model_flux{1}.biomass==0) = 0; % value = 0 if no biomass
    model_flux{1}.intl_idx(end+1) = numel(model_flux{1}.rxns);

elseif numModels == 2 % 2 Models
    % Calculate exchange flux from transport flux
    [exch_flux1,    ~    ] = trsptFlux2exchFlux(mets,flux_all_data.S,trspt_idx,flux{1}(trspt_idx,:));
    exch_flux1(abs(exch_flux1) < tol) = 0;
    [exch_flux2,exch_mets] = trsptFlux2exchFlux(mets,flux_all_data.S,trspt_idx,flux{2}(trspt_idx,:));
    [~,exch_mets_idx,~] = intersect(mets,exch_mets);
    exchRxns_metsIdx = findExchRxnsFromExchMets(temp_model,exch_mets_idx);
    exch_flux2(abs(exch_flux2) < tol) = 0;    
    % Check new exchange fluxes
    total_flux = flux{1}(exchRxns_metsIdx,:);
    calc_flux = exch_flux1+exch_flux2;
    diff_flux = total_flux - calc_flux;
    if ~isempty(find(abs(diff_flux) > tol))
        warning('myfuns:process_flux_all_new:IncorrectCalc', ...
            'Incorrect calculation of exchange flux');
    end
    flux{1}(exchRxns_metsIdx,:) = exch_flux1;
    flux{2}(exchRxns_metsIdx,:) = exch_flux2;
    
    % Sort Biomass According to Binary Variables   
    M1_biomass = biomass{1};
    M2_biomass = biomass{2};
    M1_flux = flux{1};
    M2_flux = flux{2};
    M1_int = int{1};
    M2_int = int{2};
   
    % Classify Based on Minimum Distance
    for ii = numel(M1_biomass)-1:-1:1 % start at unconstrained       
        % Calculate Distance
        D11 = pdist2(M1_int(:,ii)', M1_int(:,ii+1)'); % model 1 at this sparse_con to the previous sparse_con
        D22 = pdist2(M2_int(:,ii)', M2_int(:,ii+1)'); % model 2 at this sparse_con to the previous sparse_con
        D12 = pdist2(M1_int(:,ii)', M2_int(:,ii+1)'); % model 1 at this sparse_con to model 2 at the previous sparse_con
        D21 = pdist2(M2_int(:,ii)', M1_int(:,ii+1)'); % model 2 at this sparse_con to model 1 at the previous sparse_con
        
        % Model Prediction
        if D12 < D11 || D21 < D22 % swap
            M1_biomass(:,ii) = biomass{2}(ii);
            M2_biomass(:,ii) = biomass{1}(ii);
            M1_flux(:,ii) = flux{2}(:,ii);
            M2_flux(:,ii) = flux{1}(:,ii);
            M1_int(:,ii) = int{2}(:,ii);
            M2_int(:,ii) = int{1}(:,ii);
        else % keep the same
            M1_biomass(:,ii) = biomass{1}(ii);
            M2_biomass(:,ii) = biomass{2}(ii);
            M1_flux(:,ii) = flux{1}(:,ii);
            M2_flux(:,ii) = flux{2}(:,ii);
            M1_int(:,ii) = int{1}(:,ii);
            M2_int(:,ii) = int{2}(:,ii);
        end
    end

    % Model 1
    model_flux{1}.sparse_con = sparse_con;
    model_flux{1}.mets = mets;
    model_flux{1}.rxns = rxns;
    model_flux{1}.biomass = M1_biomass;
    model_flux{1}.flux = M1_flux;
    model_flux{1}.int = M1_int;    
    model_flux{1}.exch_idx = exch_idx;
    model_flux{1}.trspt_idx = trspt_idx;
    model_flux{1}.intl_idx = intl_idx;
    model_flux{1}.total_flux = total_flux;
    model_flux{1}.calc_flux = calc_flux;
    % Add ATPM Reaction
    model_flux{1}.rxns{end+1} = atpm_info.atpm_name;
    model_flux{1}.flux(end+1,:) = atpm_info.atpm_value;
    model_flux{1}.flux(end,model_flux{1}.biomass==0) = 0; % value = 0 if no biomass
    model_flux{1}.int(end+1,:) = 1;
    model_flux{1}.int(end,model_flux{1}.biomass==0) = 0; % value = 0 if no biomass
    model_flux{1}.intl_idx(end+1) = numel(model_flux{1}.rxns);
    
    % Model 2
    model_flux{2}.sparse_con = sparse_con;
    model_flux{2}.mets = mets;
    model_flux{2}.rxns = rxns;
    model_flux{2}.biomass = M2_biomass;
    model_flux{2}.flux = M2_flux;
    model_flux{2}.int = M2_int;  
    model_flux{2}.exch_idx = exch_idx;
    model_flux{2}.trspt_idx = trspt_idx;
    model_flux{2}.intl_idx = intl_idx;
    model_flux{2}.total_flux = total_flux;
    model_flux{2}.calc_flux = calc_flux;
    % Add ATPM Reaction
    model_flux{2}.rxns{end+1} = atpm_info.atpm_name;
    model_flux{2}.flux(end+1,:) = atpm_info.atpm_value;
    model_flux{2}.flux(end,model_flux{2}.biomass==0) = 0; % value = 0 if no biomass
    model_flux{2}.int(end+1,:) = 1;
    model_flux{2}.int(end,model_flux{2}.biomass==0) = 0; % value = 0 if no biomass
    model_flux{2}.intl_idx(end+1) = numel(model_flux{2}.rxns);
end




