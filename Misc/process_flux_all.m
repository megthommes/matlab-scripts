function [model] = process_flux_all(flux_all,sparse_con,S_info,rxns,numModels)
%%PROCESS_FLUX_ALL process flux_all into data that can be parsed and
%%created into Escher flux maps
%
% model = PROCESS_FLUX_ALL(flux_all,sparse_con,S,numModels)
%
%REQUIRED INPUT
% flux_all: output of the reaction-constraint algorithm, matrix [rxns x sparse_con]
% sparse_con: sparsity constraints, vector
% S_info: structure that must contain the following fields:
%   N_u: number of exchange/uptake reactions
%   N_e: number of transport/exchange reactions
%   N_i: number of intracellular/internal reactions
%   bio_idx: index of the biomass reaction in rxns
%   atpm_name: name of the ATPM reaction to be added
%   atpm_value: flux of the ATPM reaction
% rxns: reaction names
% numModels: number of models flux_all needs to be partitioned into
%
%OUTPUT
% model: cell structure that contains the following fields:
%   sparse_con: sparsity constraints, vector
%   rxns: reaction names, vector
%   biomass: biomass flux, vector
%   flux: reaction fluxes, matrix [rxns x sparse_con]
%   int: reaction binary variables t, matrix [rxns x sparse_con]
%
% Meghan Thommes 1/11/2017

%% Reassign Model Labels

if numModels == 1 % 1 Model
    exch_flux = 1:S_info.N_u; % exchange/uptake reaction flux
    trspt_flux = (1:S_info.N_e) + exch_flux(end); % transport/exchange reaction flux
    intl_flux = (1:S_info.N_i-1) + trspt_flux(end); % intracellular/internal reaction flux
    intl_int = (1:S_info.N_i-1) + intl_flux(end); % intracellular/internal reaction binary variable t
    
    model{1}.sparse_con = sparse_con;
    model{1}.rxns = rxns([exch_flux, trspt_flux, intl_flux]);
    model{1}.biomass = flux_all(S_info.bio_idx,:);
    model{1}.flux = flux_all([exch_flux trspt_flux intl_flux],:);
    model{1}.int = [repmat([ones(numel(exch_flux),1); ones(numel(trspt_flux),1)],1,numel(sparse_con)); flux_all(intl_int,:)];
    % Add ATPM Reaction
    model{1}.rxns{end+1} = S_info.atpm_name;
    model{1}.flux(end+1,:) = S_info.atpm_value.*ones(1,numel(sparse_con));
    model{1}.int(end+1,:) = 1;
    
    clear exch* trspt* intl*
elseif numModels == 2 % 2 Models
    exch_flux = 1:S_info.N_u; % exchange/uptake reaction flux
    trspt_flux1 = (1:S_info.N_e) + exch_flux(end); % transport/exchange reaction flux, model 1
    intl_flux1 = (1:S_info.N_i-1) + trspt_flux1(end); % intracellular/internal reaction flux, model 1
    trspt_flux2 = (1:S_info.N_e) + intl_flux1(end); % transport/exchange reaction flux, model 2
    intl_flux2 = (1:S_info.N_i-1) + trspt_flux2(end); % intracellular/internal reaction flux, model 2
    intl_int1 = (1:S_info.N_i-1) + intl_flux2(end); % intracellular/internal reaction binary variable t, model 1
    intl_int2 = (1:S_info.N_i-1) + intl_int1(end); % intracellular/internal reaction binary variable t, model 2
       
    M1_biomass = flux_all(S_info.bio_idx,:);
    M2_biomass = flux_all(S_info.N_e + S_info.N_i - 1 + S_info.bio_idx,:);
    M1_flux = flux_all([exch_flux trspt_flux1 intl_flux1],:);
    M2_flux = flux_all([exch_flux trspt_flux2 intl_flux2],:);
    M1_int = [ones(length(exch_flux),length(sparse_con)); ones(length(trspt_flux1),length(sparse_con)); flux_all(intl_int1,:)];
    M2_int = [ones(length(exch_flux),length(sparse_con)); ones(length(trspt_flux2),length(sparse_con)); flux_all(intl_int2,:)];
    
    newM1_biomass = zeros(size(M1_biomass)); newM1_biomass(:,1:2) = M1_biomass(:,1:2);
    newM2_biomass = zeros(size(M2_biomass)); newM2_biomass(:,1:2) = M2_biomass(:,1:2);
    newM1_flux = zeros(size(M1_flux)); newM1_flux(:,1) = M1_flux(:,1);
    newM2_flux = zeros(size(M2_flux)); newM2_flux(:,1) = M2_flux(:,1);
    newM1_int = zeros(size(M1_int)); newM1_int([exch_flux trspt_flux1],:) = 1; newM1_int(:,1) = M1_int(:,1);
    newM2_int = zeros(size(M2_int)); newM2_int([exch_flux trspt_flux1],:) = 1; newM2_int(:,1) = M2_int(:,1);
    
    % Calculate Euclidean Distance
    for ii = 3:size(M1_biomass,2)
        % Model 1
        D11 = pdist2(M1_biomass(:,ii-1),M1_biomass(:,ii));
        D12 = pdist2(M1_biomass(:,ii-1),M2_biomass(:,ii));
        D1 = [D11 D12];
        minD1_idx = find(D1 == min(D1),1);
        % Model 2
        D22 = pdist2(M2_biomass(:,ii-1),M2_biomass(:,ii));
        D21 = pdist2(M2_biomass(:,ii-1),M1_biomass(:,ii));
        D2 = [D22 D21];
        minD2_idx = find(D2 == min(D2),1);
        % Closest Distance
        if minD1_idx == 1 && minD2_idx == 1
            newM1_biomass(:,ii) = M1_biomass(:,ii);
            newM2_biomass(:,ii) = M2_biomass(:,ii);
            newM1_flux(:,ii) = M1_flux(:,ii);
            newM2_flux(:,ii) = M2_flux(:,ii);
            newM1_int(:,ii) = M1_int(:,ii);
            newM2_int(:,ii) = M2_int(:,ii);
        elseif minD1_idx == 2 && minD2_idx == 2
            newM1_biomass(:,ii) = M2_biomass(:,ii);
            newM2_biomass(:,ii) = M1_biomass(:,ii);
            newM1_flux(:,ii) = M2_flux(:,ii);
            newM2_flux(:,ii) = M1_flux(:,ii);
            newM1_int(:,ii) = M2_int(:,ii);
            newM2_int(:,ii) = M1_int(:,ii);
        else
            'error'
        end
        M1_biomass(:,ii) = newM1_biomass(:,ii);
        M2_biomass(:,ii) = newM2_biomass(:,ii);
        M1_flux(:,ii) = newM1_flux(:,ii);
        M2_flux(:,ii) = newM2_flux(:,ii);
        M1_int(:,ii) = newM1_int(:,ii);
        M2_int(:,ii) = newM2_int(:,ii);
    end
    % model 1
    model{1}.sparse_con = sparse_con;
    model{1}.rxns = rxns([exch_flux, trspt_flux1, intl_flux1]);
    [~,~,bio_idx] = intersect(rxns(S_info.bio_idx),model{1}.rxns);
    model{1}.biomass = M1_flux(bio_idx,:);
    model{1}.flux = M1_flux;
    model{1}.int = M1_int;
    % Add ATPM Reaction
    model{1}.rxns{end+1} = S_info.atpm_name;
    model{1}.flux(end+1,:) = S_info.atpm_value.*ones(1,numel(sparse_con));
    model{1}.int(end+1,:) = 1;
    % model 2
    model{2}.sparse_con = sparse_con;
    model{2}.rxns = rxns([exch_flux, trspt_flux1, intl_flux1]);
    [~,~,bio_idx] = intersect(rxns(S_info.bio_idx),model{2}.rxns);
    model{2}.biomass = M2_flux(bio_idx,:);
    model{2}.flux = M2_flux;
    model{2}.int = M2_int;
    % Add ATPM Reaction
    model{2}.rxns{end+1} = S_info.atpm_name;
    model{2}.flux(end+1,:) = S_info.atpm_value.*ones(1,numel(sparse_con));
    model{2}.int(end+1,:) = 1;
    clear M* newM*
    
    clear D* exch* trspt* intl*
end