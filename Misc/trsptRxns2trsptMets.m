function [trspt_mets_flux_total,trspt_mets,trspt_rxns] = trsptRxns2trsptMets(rxns,mets,S,trspt_idx,trspt_flux)
%TRSPTRXNS2TRSPTMETS Calculate the total flux of extracellular metabolites
%based on transport reaction fluxes.
%   [trspt_mets_flux_total,trspt_mets,trspt_rxns] = TRSPTRXNS2TRSPTMETS(rxns,mets,S,trspt_idx,trspt_flux)
%
% INPUTS
% rxns: Reaction names
% mets: Metabolite names
% S: Stoichiometric matrix [mets x rxns]
% trspt_idx: Indices of transport reactions
% trspt_flux: Flux of each transport reaction
%
% OUTPUTS
% trspt_flux_total: Total extracellular flux of metabolites [extracellular mets x transport rxns]
% trspt_mets: Names of extracellular metabolites (labels on trspt_flux_total rows)
% trspt_rxns: Names of transport reactions (labels on trspt_flux_total columns)
%
% Meghan Thommes 02/23/2017

exch_idx = find(full((sum(S==-1,1) == 1) & (sum(S~=0) == 1))'~=0);
% Pre-Allocate
newS = zeros(numel(exch_idx),numel(trspt_idx));
trspt_mets = cell.empty(0);
trspt_rxns = cell.empty(0);
num_mets = 0; % tracks number of extracellular metabolites

for ii = 1:numel(trspt_idx)
    trspt_rxns(ii) = rxns(trspt_idx(ii)); % transport reaction name
    mets_idx = find(full(S(:,trspt_idx(ii)))~=0); % indices of metabolites associated with the transport reaction
    mets_temp = mets(mets_idx); % names of metabolites associated with the transport reaction
    
    for jj = 1:numel(mets_idx) % for all metabolites in this transport reaction
        k = strfind(mets_temp(jj),'[e]');
        if ~isempty(k{:}) % if extracellular met
            % Check if Metabolite is in trspt_mets
            [~,idx,~] = intersect(trspt_mets,mets(mets_idx(jj)));
            if isempty(idx) % if metabolite is NOT in trspt_mets
                num_mets = num_mets + 1; % add another metabolite
                trspt_mets(num_mets) = mets(mets_idx(jj));
                newS(num_mets,ii) = full(S(mets_idx(jj),trspt_idx(ii)));
            else % if metabolite IS in trspt_mets
                trspt_mets(idx) = mets(mets_idx(jj));
                newS(idx,ii) = full(S(mets_idx(jj),trspt_idx(ii)));
            end
        end
    end
end
if numel(exch_idx) ~= numel(trspt_mets)
    disp('Length of trspt_mets is not equal to number of exchange reactions')
end

% Calculate Total Extracellular Flux of Metabolites
trspt_mets_flux_total = newS*trspt_flux;

end
