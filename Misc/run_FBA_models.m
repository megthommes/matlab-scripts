function [ output_args ] = run_FBA_models(MURI_model,MURI_flux)
%RUN_FBA_MODELS Run FBA on MURI models and compare their flux to the flux
%calculated from the algorithm
%
%INPUTS
% MURI_model: cell structure of the metabolic models for each of the
%   sparsity constraints. Contains the following fields:
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
% MURI_flux

%% FBA

minSumFlag = true;
tol = 1E-6;

color_rxns = cbrewer('div','Spectral',numel(MURI_model{1}.rxns));
marker_rxns = repmat('o',numel(MURI_model{1}.rxns),1); marker_rxns(MURI_model{1}.exch_idx) = 's'; marker_rxns(MURI_model{1}.trspt_idx) = 'd';
    
for num_sparseCon = 1:numel(MURI_model)
    sol{num_sparseCon} = FBA(MURI_model{num_sparseCon},'',minSumFlag);
    
    %% Plots
    

    
    % Find where FBA flux solution is different from the algorithm
    diff_flux = sol.fluxes - MURI_flux;
    idx = find(abs(diff_flux) > tol);
    
    
end

