function [time,biomass,mets_amt,exchange_rates,flux,ec_mets] = dFBA(model,media_names,media_mol,varargin)
%dFBA Solve a dynamic flux balance analysis problem for multiple organisms
%
% Solves LP problems of the form: max/min c'*v
%                                 subject to S*v = b
%                                            lb <= v <= ub
% [time,biomass,mets_amt,exchange_rates,flux] = dFBA(model,media_names,media_mol)
% [time,biomass,mets_amt,exchange_rates,flux] = dFBA(model,media_names,media_mol,Km,Vmax,dt,N,biomass_0,max_biomass)
%
%REQUIRED INPUTS
% The model structure must contain the following fields:
%   model.S:     Stoichiometric matrix
%   model.b:     Right hand side = dx/dt
%   model.c:     Objective coefficients
%   model.lb:    Lower bounds
%   model.ub:    Upper bounds
% media_names: cell string of extracellular media
% media_mol: vector of inital extracellular media concentrations [mmol]
%
%NOTE: I assume that the volume is 1 mL -> conc is mmol/mL
%
%OPTIONAL INPUTS - dFBA Paramters
%NOTE: Specify default using '' not ~
% Km: binding constant [mM]
%   default = 0.01
% Vmax: [mmol/gCDW/hr]
%   default = -10
% dt: time step [hr]
%   default = 0.1
% N: number of steps
%   default = 200
% biomass_0: initial biomass [gCDW]
%   default = 1E-5
% max_biomass: maximum biomass [gCDW]
%   default = 2.2E-4
%
%OUTPUT
% time: time vector [hr]
% biomass: cell vector of biomass {organism}[gCDW]
% mets_amt: matrix of extracellular media concentrations [time x mets]
% exchange_rates: cell matrix of media fluxes {organism}[time x mets]
% flux: cell matrix of all fluxes {organism}[time x rxns]
% ec_mets: vector of all the extraceullular media concentration names

% Meghan Thommes 9/3/2015

%% Check Inputs

if (nargin < 3)
    error('myfuns:dFBA:NotEnoughInputs', ...
        'Not enough inputs: need model, media_names, and media_mol');
elseif (nargin==3)
    if ~isstruct(model)
        error('myfuns:dFBA:IncorrectType', ...
            '"model" needs to be a structure');
    elseif ~isfield(model,'S') || ~isfield(model,'b') || ~isfield(model,'c') || ~isfield(model,'lb') || ~isfield(model,'ub')
        error('myfuns:dFBA:IncorrectType', ...
            '"model" needs "S", "b", "c", "lb", and "ub" fields');
    elseif ~iscellstr(media_names)
        error('myfuns:dFBA:IncorrectType', ...
            '"media_names" needs to be a string cell');
    elseif ~isnumeric(media_mol)
        error('myfuns:dFBA:IncorrectType', ...
            '"media_mol" needs to be a numeric array');
    end
end

%% Define Parameter Values

% Set defaults for optional inputs
optargs = {0.01.*ones(length(model),1) -10.*ones(length(model),1) 0.1 ...
    200 1E-5.*ones(length(model),1) 2.2E-4.*ones(length(model),1)};
% Skip new inputs if they are empty
newVals = cellfun(@(x) ~isempty(x), varargin);
% Overwrite parameters set by varargin
optargs(newVals) = varargin(newVals);
% Set values for optional inputs
[km,vmax,dt,N,biomass_0,max_biomass] = optargs{:};
clear optargs varargin
% volume
volume = 1E-3; % 1 mL

rmpath('C:\Users\mthommes\Documents\MATLAB\cobra\');

%% Define Extracellular Metabolites

ec_mets = {}; % extracellular metabolites
excRxns = cell(length(model)); % pre-allocate
excMets = cell(length(model)); % pre-allocate
for numModel = 1:length(model)
    % Find Exchange Reactions and Metabolites
    [excRxns{numModel},~] = findExchRxns(model{numModel});
    [excMets{numModel},~] = findExchMets(model{numModel});
    
    % Define All the Extracellular Metabolites
    mets = setdiff(model{numModel}.mets(excMets{numModel}),ec_mets); % model mets NOT in extracellular mets
    ec_mets = [mets; ec_mets]; % add to list
end

%% Match Exchange Reactions to Extracellular Metabolites

media = cell(size(model));
ec_rxns_idx = cell(size(model));

for numModel = 1:length(model)
    [isMet, index] = ismember(model{numModel}.mets,ec_mets);
    media{numModel} = index(isMet); % extracellular mets in model mets - index into extracellular mets
    clear isMet index
    [ec_rxns_idx{numModel},~] = findRxnsFromMets(model{numModel},ec_mets(media{numModel})); % model rxns of extracellular mets - index into model rxns
    ec_rxns_idx{numModel} = find(ec_rxns_idx{numModel} == 1);
end

%% Set Initial Metabolites

mets_amt_0 = zeros(size(ec_mets));
[~,media_names_idx,ec_mets_idx] = intersect(media_names,ec_mets,'stable');
mets_amt_0(ec_mets_idx) = media_mol(media_names_idx);

%% Initalize Vectors

time = (0:dt:N*dt)'; % time
mets_amt = zeros(N+1,length(mets_amt_0)); % metabolite abundance [mmol]
mets_amt(1,:) = mets_amt_0; % metabolite abundance at t=0 [mmol]
biomass = cell(length(model),1); % biomass [gCDW] 
exchange_rates = cell(length(model),1); % metabolite exchange flux [mmol/gCDW*hr]
Vuptake_max = cell(length(model),1); % maximum metabolite flux [mmol/gCDW*hr]
Vmax = cell(length(model),1); % Michaelis-Menten max reaction rate [mmol/gCDW/hr]
Km = cell(length(model),1); % Michaelis-Menten constant [mM]

% Set for each organism
for numModel = 1:length(model)
    biomass{numModel} = zeros(N+1,1); % biomass [gCDW]
    biomass{numModel}(1) = biomass_0(numModel); % biomass at t=0 [gCDW]
    exchange_rates{numModel} = zeros(N+1,length(mets_amt_0));  % metabolite exchange flux [mmol/gCDW*hr]
    Vuptake_max{numModel} = model{numModel}.lb(ec_rxns_idx{numModel})'; % maximum metabolite flux [mmol/gCDW*hr]
    Vmax{numModel} = vmax(numModel).*ones(1,length(mets_amt_0(media{numModel}))); % Michaelis-Menten max reaction rate [mmol/gCDW/hr]
    Km{numModel} = km(numModel).*ones(1,length(mets_amt_0(media{numModel}))); % Michaelis-Menten constant [mM]
end

%% Perform Dynamic Flux Balance Analysis

for n = 1:N
    % Randomize Order Organisms Uptake Metabolites
    eatOrder = randperm(length(model),length(model));
    
    for numModel = eatOrder
        % Calculate Uptake Rate for Exchange Reactions
        v_uptake = (Vmax{numModel}.*mets_amt(n,media{numModel}))./(Km{numModel}.*volume + mets_amt(n,media{numModel})); % v = v_max*([c]/(Km + [c]))
        
        % Check if Calculated Uptake is Greater than Max Uptake
        idx = find(v_uptake < Vuptake_max{numModel});
        v_uptake(idx) = Vuptake_max{numModel}(idx); % if greater than max, set to max uptake
        v_uptake(v_uptake > -1E-9 & v_uptake < 0) = 0;
        model{numModel}.lb(excRxns{numModel}) = v_uptake; % specify new bounds
        clear idx
        
        % Calculate Growth Rate
        FBA_solution = FBA(model{numModel});
        
        if strcmp(FBA_solution.status,'OPTIMAL') % need optimal solution
            growth_rate = FBA_solution.objectiveValue;
            exchange_rates{numModel}(n+1,media{numModel}) = FBA_solution.fluxes(excRxns{numModel});
            flux{numModel}(n+1,:) = FBA_solution.fluxes';
        else
            growth_rate = 0;
            exchange_rates{numModel}(n+1,media{numModel}) = zeros(1,sum(ec_rxns_idx{numModel}));
            flux{numModel}(n+1,:) = zeros(size(model{numModel}.rxns'));
        end
        
        % Calculate Biomass and External Metabolite Amounts
        biomass{numModel}(n+1) = (1 + growth_rate*dt).*biomass{numModel}(n); % X(t+dt) = X(t) + mu*X(t)*dt
        biomass{numModel}(biomass{numModel} > max_biomass) = max_biomass;
        biomass{numModel}(biomass{numModel} < 1E-9) = 0;
        
        % External Metabolite Amounts
        mets_amt(n+1,:) = mets_amt(n,:) + exchange_rates{numModel}(n,:).*biomass{numModel}(n).*dt; % c(t+dt) = c(t) + v*X(t)*dt
        mets_amt(mets_amt < 1E-9) = 0; % no negative concentrations
        
        clear v_uptake FBA_solution growth_rate
    end
    clear eatOrder
end

end