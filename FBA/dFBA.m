function [time,biomass,flux,exchMets_amt,exchMets_names,feasibilityFlag] = dFBA(model,media_names,media_amt,minSumFlag,biomass_0,dt,N,volume,km,vmax,max_biomass)
%dFBA Solve a dynamic flux balance analysis problem for multiple organisms
%
% Solves LP problems of the form: max/min c'*v
%                                 subject to S*v = b
%                                            lb <= v <= ub
% [time,biomass,flux,exchMets_amt,exchMets_names,feasibilityFlag] = dFBA(model,media_names,media_amt)
% [time,biomass,flux,exchMets_amt,exchMets_names,feasibilityFlag] = dFBA(model,media_names,media_amt,minSumFlag,biomass_0,dt,N,volume,km,vmax,max_biomass)
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
% minSumFlag: true if would like to minimize the sum of fluxes (default=false)
% biomass_0: initial biomass [gCDW]
%   default = 1E-5
% dt: time step [hr]
%   default = 0.1
% N: number of steps
%   default = 100
% volume: volume [L]
%   default = 1
% km: binding constant [mM]
%   default = 0.01
% vmax: uptake rate [??]
%   default = 10
% max_biomass: maximum biomass [gCDW]
%   default = 1
%
%OUTPUT
% time: time vector [hr]
% biomass: cell vector of biomass {organism}[gCDW]
% flux: cell matrix of all fluxes {organism}(time x rxns) [mmol/gCDW/hr]
% exchMets_amt: matrix of extracellular metabolite concentrations (time x mets) [mmol]
% exchMets_names: extracellular metabolite names
% feasibilityFlag: status of optimization at each time {organism}{time x 1}
%
% Meghan Thommes 05/23/2017 - Added feasibilityFlag
% Meghan Thommes 05/16/2017 - Added minSumFlag
% Meghan Thommes 05/12/2017 - v_uptake cannot be greater than the available
%                             substrate concentration
% Meghan Thommes 03/29/2017 - Renamed variables for clarity
% Meghan Thommes 03/02/2016 - Added volume as a variable
% Meghan Thommes 03/01/2016 - Updated Initializing Metabolites
% Meghan Thommes 09/03/2015

%% Check Inputs

if (nargin < 3)
    error('myfuns:dFBA:NotEnoughInputs', ...
        'Not enough inputs: need model, media_names, and media_amt');
elseif (nargin==3)
    for ii = 1:numel(model)
        if ~isstruct(model{ii})
            error('myfuns:dFBA:IncorrectType', ...
                '"model" needs to be a structure');
        elseif ~isfield(model{ii},'S') || ~isfield(model{ii},'b') || ~isfield(model{ii},'c') || ~isfield(model{ii},'lb') || ~isfield(model{ii},'ub')
            error('myfuns:dFBA:IncorrectType', ...
                '"model" needs "S", "b", "c", "lb", and "ub" fields');
        end
    end
    
    if ~iscellstr(media_names)
        error('myfuns:dFBA:IncorrectType', ...
            '"media_names" needs to be a string cell');
    elseif ~isnumeric(media_amt)
        error('myfuns:dFBA:IncorrectType', ...
            '"media_amt" needs to be a numeric array');
    end
end

%% Define Parameter Values

if ~exist('minSumFlag','var') || isempty(minSumFlag)
    minSumFlag = false; % do not minimize sum of fluxes
end
if ~exist('biomass_0','var') || isempty(biomass_0)
    biomass_0 = 1E-5.*ones(numel(model),1); % initial biomass [gCDW]
end
if ~exist('dt','var') || isempty(dt)
    dt = 0.1; % time step [hr]
end
if ~exist('N','var') || isempty(N)
    N = 100; % number of steps
end
if ~exist('volume','var') || isempty(volume)
    volume = 1; % volume [L]
end
if ~exist('km','var') || isempty(km)
    km = 0.01.*ones(numel(model),1); % binding constant [mM]
end
if ~exist('vmax','var') || isempty(vmax)
    vmax = 10.*ones(numel(model),1); % 
end
if ~exist('max_biomass','var') || isempty(max_biomass)
    max_biomass = 2.2E-4.*ones(numel(model),1); % maximum biomass [gCDW]
end

% Thresholds below which is zero
flux_thresh = 1E-9;
medium_thresh = 1E-9;
biomass_thresh = 1E-9;

%% Define Extracellular Metabolites

exchMets_names = {}; % extracellular metabolites
exchRxns_idxInModel = cell(numel(model)); % pre-allocate
exchMets_idxInModel = cell(numel(model)); % pre-allocate
for numModel = 1:numel(model)
    % Find Exchange Reactions and Metabolites
    [exchRxns_idxInModel{numModel},~] = findExchRxns(model{numModel});
    [exchMets_idxInModel{numModel},~] = findExchMets(model{numModel});
    
    % Define All the Extracellular Metabolites
    mets = setdiff(model{numModel}.mets(exchMets_idxInModel{numModel}),exchMets_names); % model mets NOT in extracellular mets
    exchMets_names = [mets; exchMets_names]; % add to list
end

%% Match Exchange Reactions to Extracellular Metabolites

exchMets_idxInMedia = cell(size(model));
exchRxns_idxInMedia = cell(size(model));

for numModel = 1:numel(model)
    [isMet, index] = ismember(model{numModel}.mets,exchMets_names);
    exchMets_idxInMedia{numModel} = index(isMet); % extracellular mets in model mets - index into extracellular mets
    clear isMet index
    [exchRxns_idxInMedia{numModel},~] = findExchRxnsFromMets(model{numModel},exchMets_names(exchMets_idxInMedia{numModel})); % model rxns of extracellular mets - index into model rxns
end

%% Set Initial Metabolites

exchMets_amt0 = zeros(size(exchMets_names));
[~,media_names_idx,exchMets_names_idx] = intersect(media_names,exchMets_names,'stable');
exchMets_amt0(exchMets_names_idx) = media_amt(media_names_idx);

%% Initalize Vectors

time = (0:dt:N*dt)'; % time
exchMets_amt = zeros(N+1,numel(exchMets_amt0)); % metabolite abundance [mmol]
exchMets_amt(1,:) = exchMets_amt0; % metabolite abundance at t=0 [mmol]
biomass = cell(numel(model),1); % biomass [gCDW] 
flux = cell(numel(model),1); % metabolite exchange flux [mmol/gCDW/hr]
Vuptake_max = cell(numel(model),1); % maximum metabolite uptake rate [mmol/gCDW/hr]
Vmax = cell(numel(model),1); % Michaelis-Menten max reaction rate [mmol/gCDW/hr]
Km = cell(numel(model),1); % Michaelis-Menten constant [mM]
feasibilityFlag = cell(numel(model),1);

% Set for each organism
for numModel = 1:numel(model)
    biomass{numModel} = zeros(N+1,1); % biomass [gCDW]
    biomass{numModel}(1) = biomass_0(numModel); % biomass at t=0 [gCDW]
    flux{numModel} = zeros(N+1,numel(model{numModel}.rxns));  % metabolite exchange flux [mmol/gCDW/hr]
    Vuptake_max{numModel} = model{numModel}.lb(exchRxns_idxInMedia{numModel})'; % maximum metabolite uptake rate [mmol/gCDW/hr]
    Vmax{numModel} = -vmax(numModel).*ones(1,numel(Vuptake_max{numModel})); % Michaelis-Menten max reaction rate [mmol/gCDW/hr]
    Km{numModel} = km(numModel).*ones(1,numel(exchMets_amt0(exchMets_idxInMedia{numModel}))); % Michaelis-Menten constant [mM]
    feasibilityFlag{numModel} = cell(N+1,1);
end

%% Perform Dynamic Flux Balance Analysis

for n = 1:N
    % Randomize Order Organisms Uptake Metabolites
    eatOrder = randperm(numel(model),numel(model));
    
    for numModel = eatOrder   
        % Calculate Uptake Rate for Exchange Reactions
        v_uptake = (Vmax{numModel}.*exchMets_amt(n,exchMets_idxInMedia{numModel}))./(Km{numModel}.*volume + exchMets_amt(n,exchMets_idxInMedia{numModel})); % v = v_max*(S(t)/(Km + S(t)))
        
        % Check if Calculated Uptake is Greater than Available Metabolite Concentration or Max Uptake Rate
        v_uptake(v_uptake < Vuptake_max{numModel}) = Vuptake_max{numModel}(v_uptake < Vuptake_max{numModel}); % if greater than max, set to max uptake
        exchMets_idxGreaterThanConc = v_uptake.*biomass{numModel}(n).*dt > exchMets_amt(n,:);
        v_uptake(exchMets_idxGreaterThanConc) = -exchMets_amt(n,exchMets_idxGreaterThanConc)./(biomass{numModel}(n).*dt); % v = -S(t)/(X(t)*dt)
        
        % Specify New Bounds
        v_uptake(v_uptake > -flux_thresh & v_uptake < 0) = 0;
        model{numModel}.lb(exchRxns_idxInModel{numModel}) = v_uptake;
        
        % Calculate Growth Rate
        FBA_solution = FBA(model{numModel},'',minSumFlag);
        
        if strcmp(FBA_solution.status,'OPTIMAL') % need optimal solution
            growth_rate = FBA_solution.objectiveValue;
            flux{numModel}(n+1,:) = FBA_solution.fluxes';
            exchMets_flux = FBA_solution.fluxes(exchRxns_idxInModel{numModel})';
            feasibilityFlag{numModel}{n+1} = FBA_solution.status;
        else
            growth_rate = 0;
            flux{numModel}(n+1,:) = zeros(1,numel(model{numModel}.rxns));
            exchMets_flux = zeros(1,numel(exchRxns_idxInMedia{numModel}));
            feasibilityFlag{numModel}{n+1} = FBA_solution.status;
        end
        
        % Calculate Biomass and External Metabolite Amounts
        biomass{numModel}(n+1) = (1 + growth_rate*dt).*biomass{numModel}(n); % X(t+dt) = X(t) + mu*X(t)*dt
        biomass{numModel}(biomass{numModel} > max_biomass) = max_biomass;
        biomass{numModel}(biomass{numModel} < biomass_thresh) = 0;
        
        % External Metabolite Amounts
        exchMets_amt(n+1,:) = exchMets_amt(n,:) + exchMets_flux.*biomass{numModel}(n).*dt; % S(t+dt) = S(t) + v*X(t)*dt
        exchMets_amt(exchMets_amt < medium_thresh) = 0; % no negative concentrations
        
        clear v_uptake FBA_solution growth_rate exchMets_flux
    end
    clear eatOrder
end

end