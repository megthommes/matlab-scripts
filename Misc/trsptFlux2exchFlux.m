function [exch_flux,exch_mets] = trsptFlux2exchFlux(mets,S,trspt_idx,trspt_flux,rxns)
%trsptFlux2exchFlux Calculate the total extracellular flux based on
%transport reaction flux.
%   [exch_flux,exch_mets] = trsptFlux2exchFlux(mets,S,trspt_idx,trspt_flux)
%
% INPUTS
% mets: Metabolite names
% S: Stoichiometric matrix [mets x rxns]
% trspt_idx: Indices of transport reactions
% trspt_flux: Flux of each transport reaction [trspt_rxns x time]
%
% OUTPUTS
% exch_flux: Extracellular flux of metabolites [exch_mets x time]
% exch_mets: Names of extracellular metabolites (labels on exch_flux rows)
%
% Meghan Thommes 07/14/2017

%% Check Inputs

if (nargin < 4)
    error('myfuns:trsptFlux2exchFlux:IncorrectInput', ...
            'Not enough inputs');
end

if numel(mets) ~= size(S,1)
    error('myfuns:trsptFlux2exchFlux:IncorrectInput', ...
        'Number of mets must be the same size as the 1st dim of S');
end

if size(trspt_flux,1) ~= numel(trspt_idx)
    trspt_flux = trspt_flux';
    if size(trspt_flux,1) ~= numel(trspt_idx)
        error('myfuns:trsptFlux2exchFlux:IncorrectInput', ...
            'Number of trspt_idx must be the same size as the 1st dim of trspt_flux');
    end
end

%% Calculate Exchange Flux

model.S = sparse(S);
[exchMets_idx,~] = findExchMets(model); % find exchange metabolite indices
exch_mets = mets(exchMets_idx);

% Pre-Allocate
newS = zeros(numel(exchMets_idx),numel(trspt_idx));

for ii = 1:numel(trspt_idx) % for each transport reaction
    mets_idx = find(S(:,trspt_idx(ii))); % indices of metabolite(s) associated with the transport reaction
    mets_idx = intersect(mets_idx,exchMets_idx); % check that the metabolite is an extracellular metabolite
    [~,idx,~] = intersect(exch_mets,mets(mets_idx)); % find index/indices in exch_mets
    newS(idx,ii) = full(S(mets_idx,trspt_idx(ii))); % update stoichimetric matrix
    clear mets_idx idx
end

% Calculate Total Extracellular Flux of Metabolites
exch_flux = newS*trspt_flux;

end

