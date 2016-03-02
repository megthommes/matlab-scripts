function [time,biomass,exchMets_amt,exchMets_flux,exchMets_names] = dFBA(model,media_names,media_amt,varargin)
%dFBA Solve a dynamic flux balance analysis problem for multiple organisms
%
% Solves LP problems of the form: max/min c'*v
%                                 subject to S*v = b
%                                            lb <= v <= ub
% [time,biomass,exchMets_amt,exchMets_flux,flux,exchMets_names] = dFBA(model,media_names,media_amt)
% [time,biomass,exchMets_amt,exchMets_flux,flux,exchMets_names] = dFBA(model,media_names,media_amt,biomass_0,dt,N,volume,Km,Vmax,max_biomass)
%
%REQUIRED INPUTS
% The model structure must contain the following fields:
%   model.S:     Stoichiometric matrix
%   model.b:     Right hand side = dx/dt
%   model.c:     Objective coefficients
%   model.lb:    Lower bounds
%   model.ub:    Upper bounds
% media_names: cell string of extracellular media
% media_amt: vector of inital extracellular media concentrations [mmol]
%
%OPTIONAL INPUTS - dFBA Paramters
%NOTE: Specify default using '' not ~
% biomass_0: initial biomass [gCDW]
%   default = 1E-5
% dt: time step [hr]
%   default = 0.1
% N: number of steps
%   default = 100
% volume: volume [L]
%   default = 1
% Km: binding constant [mM]
%   default = 0.01
% Vmax: [mmol/gCDW/hr]
%   default = -10
% max_biomass: maximum biomass [gCDW]
%   default = 1
%
%OUTPUT
% time: time vector [hr]
% biomass: cell vector of biomass {organism}[gCDW]
% exchMets_amt: matrix of extracellular media concentrations (time x mets) [mmol]
% exchMets_flux: cell matrix of media fluxes {organism}(time x mets) [mmol/gCDW/hr]
% flux: cell matrix of all fluxes {organism}(time x rxns) [mmol/gCDW/hr]
% exchMets_names: vector of all the extraceullular media concentration names
%
% Meghan Thommes 03/02/2016 - Added volume as a variable
% Meghan Thommes 03/01/2016 - Updated Initializing Metabolites
% Meghan Thommes 09/03/2015

%% Check Inputs

if (nargin < 3)
    error('myfuns:dFBA:NotEnoughInputs', ...
        'Not enough inputs: need model, media_names, and media_amt');
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
    elseif ~isnumeric(media_amt)
        error('myfuns:dFBA:IncorrectType', ...
            '"media_amt" needs to be a numeric array');
    end
end

%% Define Parameter Values

% Set defaults for optional inputs
optargs = {1E-5.*ones(length(model),1), 0.1, 100, 1, ... % biomass_0, dt, N, volume
    0.01.*ones(length(model),1), -10.*ones(length(model),1), 2.2E-4.*ones(length(model),1)}; %  Km, Vmax, max_biomass
% Skip new inputs if they are empty
newVals = cellfun(@(x) ~isempty(x), varargin);
% Overwrite parameters set by varargin
optargs(newVals) = varargin(newVals);
% Set values for optional inputs
[biomass_0,dt,N,volume,km,vmax,max_biomass] = optargs{:};
clear optargs varargin
% Thresholds below which is zero
media_thresh = 1E-9;
biomass_thresh = 1E-9;

%% Define Extracellular Metabolites

exchMets_names = {}; % extracellular metabolites
excRxns = cell(length(model)); % pre-allocate
excMets = cell(length(model)); % pre-allocate
for numModel = 1:length(model)
    % Find Exchange Reactions and Metabolites
    [excRxns{numModel},~] = findExchRxns(model{numModel});
    [excMets{numModel},~] = findExchMets(model{numModel});
    
    % Define All the Extracellular Metabolites
    mets = setdiff(model{numModel}.mets(excMets{numModel}),exchMets_names); % model mets NOT in extracellular mets
    exchMets_names = [mets; exchMets_names]; % add to list
end

%% Match Exchange Reactions to Extracellular Metabolites

exchMets_idx = cell(size(model));
ec_rxns_idx = cell(size(model));

for numModel = 1:length(model)
    [isMet, index] = ismember(model{numModel}.mets,exchMets_names);
    exchMets_idx{numModel} = index(isMet); % extracellular mets in model mets - index into extracellular mets
    clear isMet index
    [ec_rxns_idx{numModel},~] = findRxnsFromMets(model{numModel},exchMets_names(exchMets_idx{numModel})); % model rxns of extracellular mets - index into model rxns
    ec_rxns_idx{numModel} = find(ec_rxns_idx{numModel} == 1);
end

%% Set Initial Metabolites

exchMets_amt_0 = zeros(size(exchMets_names));
[~,media_names_idx,exchMets_names_idx] = intersect(media_names,exchMets_names,'stable');
exchMets_amt_0(exchMets_names_idx) = media_amt(media_names_idx);

%% Initalize Vectors

time = (0:dt:N*dt)'; % time
exchMets_amt = zeros(N+1,length(exchMets_amt_0)); % metabolite abundance [mmol]
exchMets_amt(1,:) = exchMets_amt_0; % metabolite abundance at t=0 [mmol]
biomass = cell(length(model),1); % biomass [gCDW] 
exchMets_flux = cell(length(model),1); % metabolite exchange flux [mmol/gCDW*hr]
Vuptake_max = cell(length(model),1); % maximum metabolite flux [mmol/gCDW*hr]
Vmax = cell(length(model),1); % Michaelis-Menten max reaction rate [mmol/gCDW/hr]
Km = cell(length(model),1); % Michaelis-Menten constant [mM]

% Set for each organism
for numModel = 1:length(model)
    biomass{numModel} = zeros(N+1,1); % biomass [gCDW]
    biomass{numModel}(1) = biomass_0(numModel); % biomass at t=0 [gCDW]
    exchMets_flux{numModel} = zeros(N+1,length(exchMets_amt_0));  % metabolite exchange flux [mmol/gCDW*hr]
    Vuptake_max{numModel} = model{numModel}.lb(ec_rxns_idx{numModel})'; % maximum metabolite flux [mmol/gCDW*hr]
    Vmax{numModel} = vmax(numModel).*ones(1,length(exchMets_amt_0(exchMets_idx{numModel}))); % Michaelis-Menten max reaction rate [mmol/gCDW/hr]
    Km{numModel} = km(numModel).*ones(1,length(exchMets_amt_0(exchMets_idx{numModel}))); % Michaelis-Menten constant [mM]
end

%% Perform Dynamic Flux Balance Analysis

for n = 1:N
    % Randomize Order Organisms Uptake Metabolites
    eatOrder = randperm(length(model),length(model));
    
    for numModel = eatOrder   
        % Calculate Uptake Rate for Exchange Reactions
        v_uptake = (Vmax{numModel}.*exchMets_amt(n,exchMets_idx{numModel}))./(Km{numModel}.*volume + exchMets_amt(n,exchMets_idx{numModel})); % v = v_max*([S]/(Km + [S]))
        
        % Check if Calculated Uptake is Greater than Max Uptake
        v_uptake(v_uptake < Vuptake_max{numModel}) = Vuptake_max{numModel}(v_uptake < Vuptake_max{numModel}); % if greater than max, set to max uptake
        v_uptake(v_uptake > -1E-9 & v_uptake < 0) = 0;
        model{numModel}.lb(excRxns{numModel}) = v_uptake; % specify new bounds
        
        % Calculate Growth Rate
        FBA_solution = FBA(model{numModel});
        
        if strcmp(FBA_solution.status,'OPTIMAL') % need optimal solution
            growth_rate = FBA_solution.objectiveValue;
            exchMets_flux{numModel}(n+1,exchMets_idx{numModel}) = FBA_solution.fluxes(excRxns{numModel});
        else
            growth_rate = 0;
            exchMets_flux{numModel}(n+1,exchMets_idx{numModel}) = zeros(1,length(ec_rxns_idx{numModel}));
        end
        
        % Calculate Biomass and External Metabolite Amounts
        biomass{numModel}(n+1) = (1 + growth_rate*dt).*biomass{numModel}(n); % X(t+dt) = X(t) + mu*X(t)*dt
        biomass{numModel}(biomass{numModel} > max_biomass) = max_biomass;
        biomass{numModel}(biomass{numModel} < biomass_thresh) = 0;
        
        % External Metabolite Amounts
        exchMets_amt(n+1,:) = exchMets_amt(n,:) + exchMets_flux{numModel}(n,:).*biomass{numModel}(n).*dt; % S(t+dt) = S(t) + v*X(t)*dt
        exchMets_amt(exchMets_amt < media_thresh) = 0; % no negative concentrations
        
        clear v_uptake FBA_solution growth_rate
    end
    clear eatOrder
end

end