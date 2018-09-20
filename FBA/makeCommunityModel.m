function [communityModel] = makeCommunityModel(individualModels,equalFlag)
%makeCommunityModel Make a community metabolic model from individual models
%
% The stoichiometric matrix becomes the form:
%           S_e S_te1      ... S_tej
%               S_ti1 S_i1     
%                          ... S_tij S_ij
% where
%   S_e : exchange reactions and extracellular metabolites
%   S_te: transport reactions and extracellular metabolites
%   S_ti: transport reactions and intracellular metabolites
%   S_i : intracellular reactions and intracellular metabolites
%
%REQUIRED INPUT
% individualModels is a cell structure of model structuress. Each model
% structure must contain the following fields:
%   model.S:       Stoichiometric matrix
%   model.c:       Objective - must be biomass function
%   model.lb:      Lower bounds
%   model.ub:      Upper bounds
%
%OPTIONAL INPUT
% equalFlag: true if would like to constrain each model to have the same
% growth rate (default=false)
%
%OUTPUT
% The communityModel structure will contain the following fields:
%   model.S:       Stoichiometric matrix
%   model.b:       Right hand side = dx/dt
%   model.c:       Objective
%   model.lb:      Lower bounds
%   model.ub:      Upper bounds
%   model.orgMets: 0 if all organisms, otherwise organism number as defined
%                by individualModels
%   model.orgRxns: 0 if all organisms, otherwise organism number as defined
%                by individualModels
%
% Meghan Thommes 09/20/2018

%% Check Inputs

if (nargin < 1)
    error('myfuns:FBA:NotEnoughInputs', ...
        'Not enough inputs: need individualModels');
end

numModels = numel(individualModels);
if numModels < 2
    error('myfuns:FBA:IncorrectType', ...
        'individualModels needs to contain at least 2 models');
end

for ii = 1:numModels
    if ~isstruct(individualModels{ii})
        error('myfuns:FBA:IncorrectType', ...
            ['"individualModels' int2str(ii) '" needs to be a structure']);
    elseif ~isfield(individualModels{ii},'S') || ~isfield(individualModels{ii},'c') || ~isfield(individualModels{ii},'lb') || ~isfield(individualModels{ii},'ub')
        error('myfuns:FBA:IncorrectType', ...
            ['"individualModels' int2str(ii) 'needs "S", "c", "lb", and "ub" fields']);
    end
end

if ~exist('equalFlag','var') || isempty(equalFlag)
    equalFlag = false;
end

%% Reaction and Metabolite Indices

% Pre-Allocate
% Reaction Indices
exchRxns_idx = cell(numModels,1); % exchange reaction indices
trsptRxns_idx = cell(numModels,1); % transport reaction indices
intlRxns_idx = cell(numModels,1); % intracellular reaction indices
% Metabolite Indices
exchMets_idx = cell(numModels,1); % extracellular metabolite indices
intlMets_idx = cell(numModels,1); % intracellular metabolite indices
% Number of Reactions
N_t = nan(numModels,1); % transport reactions
N_i = nan(numModels,1); % intracellular reactions
% Number of Metabolites
M_i = nan(numModels,1); % intracellular metabolites
% Biomass Reaction Indices
bio_idx = nan(numModels,1); % biomass reaction indices
% Exchange Reaction Names
exchRxns = [];
exchRxnNames = [];
% Extracellular Metabolite Names
exchMets = [];
exchMetNames = [];

for ii = 1:numModels
    % Reaction Indices
    [exchRxns_idx{ii},~] = findExchRxns(individualModels{ii});
    trsptRxns_idx{ii} = findTrsptRxns(individualModels{ii});
    intlRxns_idx{ii} = (1:numel(individualModels{ii}.rxns))'; intlRxns_idx{ii}([exchRxns_idx{ii}; trsptRxns_idx{ii}]) = [];
    % Metabolite Indices
    [exchMets_idx{ii},~] = findExchMets(individualModels{ii});
    intlMets_idx{ii} = (1:numel(individualModels{ii}.mets))'; intlMets_idx{ii}(exchMets_idx{ii}) = [];
    % Number of Reactions
    N_t(ii) = numel(trsptRxns_idx{ii});
    N_i(ii) = numel(intlRxns_idx{ii});
    % Number of Metabolites
    M_i(ii) = numel(intlMets_idx{ii});
    % Biomass Reaction Indices
    bio_idx(ii) = find(individualModels{ii}.c);
end

% Number of Reactions & Metabolites cont'd
for ii = 1:numModels
    exchRxns = [exchRxns; individualModels{ii}.rxns(exchRxns_idx{ii})];
    exchRxnNames = [exchRxnNames; individualModels{ii}.rxnNames(exchRxns_idx{ii})];
    exchMets = [exchMets; individualModels{ii}.mets(exchMets_idx{ii})];
    exchMetNames = [exchMetNames; individualModels{ii}.metNames(exchMets_idx{ii})];
end
[exchRxns,idx,~] = unique(exchRxns,'stable');
N_e = numel(exchRxns);
exchRxnNames = exchRxnNames(idx);
[exchMets,idx,~] = unique(exchMets,'stable');
M_e = numel(exchMets);
exchMetNames = exchMetNames(idx);

%% Create S, lb, ub, c, mets, rxns

% Pre-Allocate
mets = cell(M_e+sum(M_i),1);
metNames = cell(M_e+sum(M_i),1);
rxns = cell(N_e+sum(N_t)+sum(N_i),1);
rxnNames = cell(N_e+sum(N_t)+sum(N_i),1);
S = zeros(numel(mets), numel(rxns));
lb = nan(size(rxns));
ub = nan(size(rxns));
c = nan(size(rxns));
orgMets = nan(size(mets));
orgRxns = nan(size(rxns));

% S_e
S_e_ind = cell(numModels,1);
lb_e_ind = zeros(N_e,numModels);
ub_e_ind = zeros(N_e,numModels);
for ii = 1:numModels
    [~,mets_idx,~] = intersect(exchMets,individualModels{ii}.mets(exchMets_idx{ii}),'stable');
    [~,rxns_idx,~] = intersect(exchRxns,individualModels{ii}.rxns(exchRxns_idx{ii}),'stable');
    S_e_ind{ii} = zeros(numel(exchMets),numel(exchRxns));
    S_e_ind{ii}(mets_idx,rxns_idx) = individualModels{ii}.S(exchMets_idx{ii},exchRxns_idx{ii});
    lb_e_ind(rxns_idx,ii) = individualModels{ii}.lb(exchRxns_idx{ii});
    ub_e_ind(rxns_idx,ii) = individualModels{ii}.ub(exchRxns_idx{ii});
end
% Check for Differences/Unique Extracellular Reactions
S_e = S_e_ind{1};
C = combnk(1:numModels,2);
for ii = 2:numModels
    [~,r] = find(S_e ~= S_e_ind{ii});
    S_e(:,r) = S_e_ind{ii}(:,r);
end
% lb, ub
lb_e = min(lb_e_ind,[],2);
ub_e = max(ub_e_ind,[],2);
% Update S, lb, ub, c, mets, rxns
mets_idx = 1:M_e;
rxns_idx = 1:N_e;
S(mets_idx,rxns_idx) = S_e;
c(rxns_idx) = 0;
lb(rxns_idx) = lb_e;
ub(rxns_idx) = ub_e;
mets(mets_idx) = exchMets;
metNames(mets_idx) = exchMetNames;
rxns(rxns_idx) = exchRxns;
rxnNames(rxns_idx) = exchRxnNames;
orgMets(mets_idx) = 0;
orgRxns(rxns_idx) = 0;

% S_te
mets_idx = 1:M_e;
ii = 1;
rxns_idx = N_e + (1:N_t(ii));
S(mets_idx,rxns_idx) = individualModels{ii}.S(exchMets_idx{ii},trsptRxns_idx{ii});
lb(rxns_idx) = individualModels{ii}.lb(trsptRxns_idx{ii});
ub(rxns_idx) = individualModels{ii}.ub(trsptRxns_idx{ii});
c(rxns_idx) = individualModels{ii}.c(trsptRxns_idx{ii});
for ii = 2:numModels
    % S_t
    rxns_idx = N_e + sum(N_t(1:ii-1)) + sum(N_i(1:ii-1)) + (1:N_t(ii));
    S(mets_idx,rxns_idx) = individualModels{ii}.S(exchMets_idx{ii},trsptRxns_idx{ii});    
    lb(rxns_idx) = individualModels{ii}.lb(trsptRxns_idx{ii});
    ub(rxns_idx) = individualModels{ii}.ub(trsptRxns_idx{ii});
    c(rxns_idx) = individualModels{ii}.c(trsptRxns_idx{ii});
end
% S_ti
ii = 1;
mets_idx = M_e + (1:M_i);
rxns_idx = N_e + (1:N_t(ii));
S(mets_idx,rxns_idx) = individualModels{ii}.S(intlMets_idx{ii},trsptRxns_idx{ii});
for ii = 2:numModels
    % S_t
    mets_idx = M_e + sum(M_i(1:ii-1)) + (1:M_i);
    rxns_idx = N_e + sum(N_t(1:ii-1)) + sum(N_i(1:ii-1)) + (1:N_t(ii));
    S(mets_idx,rxns_idx) = individualModels{ii}.S(intlMets_idx{ii},trsptRxns_idx{ii});    
end

% S_i
ii = 1;
mets_idx = M_e + (1:M_i);
rxns_idx = N_e + N_t(1) + (1:N_i);
S(mets_idx,rxns_idx) = individualModels{ii}.S(intlMets_idx{ii},intlRxns_idx{ii});
lb(rxns_idx) = individualModels{ii}.lb(intlRxns_idx{ii});
ub(rxns_idx) = individualModels{ii}.ub(intlRxns_idx{ii});
c(rxns_idx) = individualModels{ii}.c(intlRxns_idx{ii});
mets(mets_idx) = individualModels{ii}.mets(intlMets_idx{ii});
if ~isempty(intersect(fieldnames(individualModels{ii}),'metNames'))
    metNames(mets_idx) = individualModels{ii}.metNames(intlMets_idx{ii});
end
rxns(rxns_idx) = individualModels{ii}.rxns(intlRxns_idx{ii});
if ~isempty(intersect(fieldnames(individualModels{ii}),'rxnNames'))
    rxnNames(rxns_idx) = individualModels{ii}.rxnNames(intlRxns_idx{ii});
end
orgMets(mets_idx) = ii;
orgRxns(rxns_idx) = ii;
for ii = 2:numModels
    % S_i
    mets_idx = M_e + sum(M_i(1:ii-1)) + (1:M_i);
    rxns_idx = N_e + sum(N_t(1:ii)) + sum(N_i(1:ii-1)) + (1:N_i);
    S(mets_idx,rxns_idx) = individualModels{ii}.S(intlMets_idx{ii},intlRxns_idx{ii});
    lb(rxns_idx) = individualModels{ii}.lb(intlRxns_idx{ii});
    ub(rxns_idx) = individualModels{ii}.ub(intlRxns_idx{ii});
    c(rxns_idx) = individualModels{ii}.c(intlRxns_idx{ii});
    mets(mets_idx) = individualModels{ii}.mets(intlMets_idx{ii});
    if ~isempty(intersect(fieldnames(individualModels{ii}),'metNames'))
        metNames(mets_idx) = individualModels{ii}.metNames(intlMets_idx{ii});
    end
    rxns(rxns_idx) = individualModels{ii}.rxns(intlRxns_idx{ii});
    if ~isempty(intersect(fieldnames(individualModels{ii}),'rxnNames'))
        rxnNames(rxns_idx) = individualModels{ii}.rxnNames(intlRxns_idx{ii});
    end
    orgMets(mets_idx) = ii;
    orgRxns(rxns_idx) = ii;
end

% Add Equality Constraints
if equalFlag == 1
    biomass_id = find(c);
    for ii = 2:numModels
        S(end+1,biomass_id([1,ii])) = [1,-1];
    end
end

%% Make Community Model

communityModel.S = sparse(S);
communityModel.b = zeros(size(S,1),1);
communityModel.c = c;
communityModel.lb = lb;
communityModel.ub = ub;
communityModel.mets = mets;
communityModel.metNames = metNames;
communityModel.rxns = rxns;
communityModel.rxnNames = rxnNames;
communityModel.orgMets = orgMets;
communityModel.orgRxns = orgRxns;

end

