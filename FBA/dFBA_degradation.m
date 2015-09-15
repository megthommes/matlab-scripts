function [time,biomass,mets_amt,exchange_rates,flux,ec_mets] = dFBA_degradation(model,media_names,media_mol,varargin)
%dFBA_degradation Solve a dynamic flux balance analysis problem for multiple organisms
%   and performs extracellular cellulose degradation
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
volume = 1; %1E-3; % 1 mL
% model that produces the extacellular enzymes
enzymeModel = 1;

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

% Add Cellulose & Enzyme
ec_mets{end+1} = 'cellulose[e]';
cellulose_met_idx = length(ec_mets);
ec_mets{end+1} = 'enzyme[e]';
enzyme_met_idx = length(ec_mets);

%% Match Exchange Reactions to Extracellular Metabolites

for numModel = 1:length(model)
    [isMet, index] = ismember(model{numModel}.mets,ec_mets);
    media{numModel} = index(isMet); % extracellular mets in model mets - index into extracellular mets
    clear isMet index
    [ec_rxns_idx{numModel},~] = findRxnsFromMets(model{numModel},ec_mets(media{numModel})); % model rxns of extracellular mets - index into model rxns
    ec_rxns_idx{numModel} = find(ec_rxns_idx{numModel} == 1);
end

%% Cellulose Degradation 

% Parameters
kcat_glc = 2.5E3; % turnover rate [1/hr]
Km_glc = 1.2; % [mM]
Km_enzyme = 1; % [mM]
Y_glc = 10; % glucose yield coefficient [mmol glucose/mmol cellulose]
beta = 10; % [mmol enzyme/(gCDW*hr)]
alpha = 0; % enzyme production cost [gCDW/mmol enzyme]

% Cellulose Initial Amount
cellulose_amt_0 = 1; % [mmol]

% Glucose
[isGlc, index] = ismember('glc-D[e]',ec_mets);
glc_met_idx = index(isGlc);
clear isGlc index

% Water
[isH2O, index] = ismember('h2o[e]',ec_mets);
h2o_met_idx = index(isH2O);
clear isH2O index

%% Set Initial Metabolites

mets_amt_0 = zeros(size(ec_mets));
[isMet, index] = ismember(media_names,ec_mets); % media in extracellular mets
mets_amt_0(index(isMet)) = media_mol;
clear isMet index

%% Initalize Vectors

% dFBA
time = [0:dt:N*dt]'; % time
mets_amt = zeros(N+1,length(mets_amt_0)); % metabolite abundance [mmol]
mets_amt(1,:) = mets_amt_0; % metabolite abundance at t=0 [mmol]
biomass = cell(length(model),1); % biomass [gCDW] 
exchange_rates = cell(length(model),1); % metabolite exchange flux [mmol/gCDW*hr]
Vuptake_max = cell(length(model),1); % maximum metabolite flux [mmol/gCDW*hr]
Vmax = cell(length(model),1); % Michaelis-Menten max reaction rate [mmol/gCDW/hr]
Km = cell(length(model),1); % Michaelis-Menten constant [mM]

% dFBA - Set for each organism
for numModel = 1:length(model)
    biomass{numModel} = zeros(N+1,1); % biomass [gCDW]
    biomass{numModel}(1) = biomass_0(numModel); % biomass at t=0 [gCDW]
    exchange_rates{numModel} = zeros(N+1,length(mets_amt_0));  % metabolite exchange flux [mmol/gCDW*hr]
    Vuptake_max{numModel} = model{numModel}.lb(ec_rxns_idx{numModel})'; % maximum metabolite flux [mmol/gCDW*hr]
    Vmax{numModel} = vmax(numModel).*ones(1,length(mets_amt_0(media{numModel}))); % Michaelis-Menten max reaction rate [mmol/gCDW/hr]
    Km{numModel} = km(numModel).*ones(1,length(mets_amt_0(media{numModel}))); % Michaelis-Menten constant [mM]
end

% Cellulose
mets_amt(1,cellulose_met_idx) = cellulose_amt_0; % cellulose abundance at t=0 [mmol]
mets_amt(1,glc_met_idx) = 1E-3; % glucose abundance at t=0 [mmol]

%% Perform Dynamic Flux Balance Analysis

glc_amt = mets_amt(1,glc_met_idx);
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
        model{numModel}.lb(ec_rxns_idx{numModel}) = v_uptake; % specify new bounds
        clear idx
        
        % Calculate Growth Rate
        FBA_solution = FBA(model{numModel});
        
        if strcmp(FBA_solution.status,'OPTIMAL') % need optimal solution
            growth_rate = FBA_solution.objectiveValue;
            exchange_rates{numModel}(n+1,media{numModel}) = FBA_solution.fluxes(ec_rxns_idx{numModel});
            flux{numModel}(n+1,:) = FBA_solution.fluxes';
        else
            growth_rate = 0;
            exchange_rates{numModel}(n+1,media{numModel}) = zeros(1,length(ec_rxns_idx{numModel}));
            flux{numModel}(n+1,:) = zeros(size(model{numModel}.rxns'));
        end
        
        % Calculate External Metabolite Amounts
        mets_amt(n+1,:) = mets_amt(n,:) + exchange_rates{numModel}(n,:).*biomass{numModel}(n).*dt; % c(t+dt) = c(t) + v*X(t)*dt
        
        % Calculate Biomass
        if numModel == enzymeModel
            biomass{numModel}(n+1) = (1 + growth_rate*dt).*biomass{numModel}(n) ...
                - alpha*mets_amt(n,enzyme_met_idx); % X(t+dt) = X(t) + mu*X(t)*dt - alpha*E(t)*dt
            biomass{numModel}(biomass{numModel} > max_biomass) = max_biomass;
            biomass{numModel}(biomass{numModel} < 1E-9) = 0;
            
            % Calculate Enzyme Amount
            % E(t+dt) = E(t) + b*X(t)*(C(t)/(Km_E + C(t)))
            mets_amt(n+1,enzyme_met_idx) = mets_amt(n,enzyme_met_idx) + ...
                beta*biomass{numModel}(n)*...
                (mets_amt(n,cellulose_met_idx)/(Km_enzyme + mets_amt(n,cellulose_met_idx)))*dt;
        else
            biomass{numModel}(n+1) = (1 + growth_rate*dt).*biomass{numModel}(n); % X(t+dt) = X(t) + mu*X(t)*dt
            biomass{numModel}(biomass{numModel} > max_biomass) = max_biomass;
            biomass{numModel}(biomass{numModel} < 1E-9) = 0;
        end
        
        % Calculate Glucose Amount
        % G(t+dt) = G(t) + kcat*E(t)*(C(t)/Km + C(t))*dt + v_uptake_glc*dt
        delta_glc = kcat_glc*mets_amt(n,enzyme_met_idx)*...
            (mets_amt(n,cellulose_met_idx)*mets_amt(n,h2o_met_idx))/...
            (Km_glc + mets_amt(n,cellulose_met_idx)*mets_amt(n,h2o_met_idx))*dt;
        mets_amt(n+1,glc_met_idx) = mets_amt(n+1,glc_met_idx) + delta_glc;
        
        % Calculate Cellulose Amount
        % C(t+dt) = C(t) - (1/Y_glc)*dG
        mets_amt(n+1,cellulose_met_idx) = mets_amt(n,cellulose_met_idx) + ...
            -(1/Y_glc)*delta_glc/dt;
        
        % Calculate Water Amount
        % W(t+dt) = W(t) - dG + v_uptake_glc*dt
        mets_amt(n+1,h2o_met_idx) = mets_amt(n+1,h2o_met_idx) - delta_glc/dt;
        
        % No Negative Abundances
        mets_amt(mets_amt < 1E-9) = 0; 
        
        clear v_uptake FBA_solution growth_rate delta_glc
    end
    clear eatOrder
    
    % Total Glucose Produced
    glc_amt = glc_amt + mets_amt(n+1,glc_met_idx);
end

end