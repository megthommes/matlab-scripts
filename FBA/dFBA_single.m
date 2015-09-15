function [time,biomass,mets_conc,exchange_rates] = dFBA(model,media_names,media_mol,varargin)
%dFBA Solve a dynamic flux balance analysis problem
%
% Solves LP problems of the form: max/min c'*v
%                                 subject to S*v = b
%                                            lb <= v <= ub
% [time,biomass,mets_conc,exchange_rates] = dFBA(model,media_names,media_mol)
% [time,biomass,mets_conc,exchange_rates] = dFBA(model,media_names,media_mol,Km,Vmax,dt,N,biomass_0,max_biomass)
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
% Vmax: [mmol/g/hr]
%   default = -10
% dt: time step [hours]
%   default = 0.1
% N: number of steps
%   default = 200
% biomass_0: initial biomass [g]
%   default = 1E-5
% max_biomass: maximum biomass [g]
%   default = 2.2E-4
%
%OUTPUT
% time: time vector [hr]
% biomass: vector of biomass [gCDW]
% mets_conc: matrix of extracellular media concentrations [mets x time]
% exchange_rates: matrix of media fluxes [mets x time]
%
% Meghan Thommes 6/11/2015 - No longer exit when infeasible, just set all
%   fluxes to zero
% Meghan Thommes 5/26/2015 - Removed conversion factor, updated setting
%   initial metabolites, made Vmax negative, & changed maximum uptake rate
%   to initial lower bound value instead of Vmax
% Meghan Thommes 4/15/2015 - Added conversion factor to exchange_rates so
%   that it has units [mmol/(gDW*hr)]
% Meghan Thommes 3/20/2015
% Meghan Thommes 3/3/2015

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
optargs = {0.01 -10 0.1 200 1E-5 2.2E-4};
% Skip new inputs if they are empty
newVals = cellfun(@(x) ~isempty(x), varargin);
% Overwrite parameters set by varargin
optargs(newVals) = varargin(newVals);
% Set values for optional inputs
[Km,Vmax,dt,N,biomass_0,max_biomass] = optargs{:};
clear optargs varargin

if (nargin > 3)
    if (nargin >= 4)
        if ~isnumeric(Km) || ~isscalar(Km)
            error('myfuns:dFBA:IncorrectType', ...
                '"Km" needs to be a numeric scalar')
        end
        if (nargin >= 5)
            if ~isnumeric(Vmax) || ~isscalar(Vmax)
                error('myfuns:dFBA:IncorrectType', ...
                    '"Vmax" needs to be a numeric scalar')
            end
            if (nargin >= 6)
                if ~isnumeric(dt) || ~isscalar(dt)
                    error('myfuns:dFBA:IncorrectType', ...
                        '"dt" needs to be a numeric scalar')
                end
                if (nargin >= 7)
                    if ~isnumeric(N) || ~isscalar(N)
                        error('myfuns:dFBA:IncorrectType', ...
                            '"N" needs to be a numeric scalar')
                    end
                    if (nargin >= 8)
                        if ~isnumeric(biomass_0) || ~isscalar(biomass_0)
                            error('myfuns:dFBA:IncorrectType', ...
                                '"biomass_0" needs to be a numeric scalar')
                        end
                        if (nargin >= 9)
                            if ~isnumeric(max_biomass) || ~isscalar(max_biomass)
                                error('myfuns:dFBA:IncorrectType', ...
                                    '"max_biomass" needs to be a numeric scalar')
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Find Exchange Reactions and Metabolites

% Exchange Reactions
[excRxns,~] = findExchRxns(model);

%% Set Initial Metabolites

[rxns_idx,~] = findRxnsFromMets(model,media_names); % index in model.rxns
[isRxn, index] = ismember(find(rxns_idx~=0),find(excRxns~=0));
excRxns = find(excRxns~=0);
idx = excRxns(index(isRxn))'; % index in excRxns
idx = idx - (max(idx) - length(excRxns));
mets_conc_0 = zeros(size(excRxns));
mets_conc_0(idx) = media_mol;
clear idx

% Check and make sure that the metabolite names match those given

%% Perform Dynamic Flux Balance Analysis - Step 1

% Initialize Vectors
Vuptake_max = model.lb(excRxns);
time = [0:dt:N*dt]'; % time
biomass = zeros(N,length(biomass_0)); % biomass
mets_conc = zeros(N,length(mets_conc_0)); % concentration
exchange_rates = zeros(N,length(mets_conc_0)); % concentration

% Calculate Flux for Exchange Reactions
v_uptake = (Vmax.*mets_conc_0)./(Km + mets_conc_0); % v = v_max*([c]/(Km + [c]))

% Check if calculated uptake is greater than max
idx = find(v_uptake < Vuptake_max);
v_uptake(idx) = Vuptake_max(idx); % if greater than max, set to max uptake
v_uptake(v_uptake > -1E-9 & v_uptake < 0) = 0;
model.lb(excRxns) = v_uptake; % specify new bounds

% FBA
FBA_solution = FBA(model);

if strcmp(FBA_solution.status,'OPTIMAL') % need optimal solution
    growth_rate = FBA_solution.objectiveValue;
    exchange_rates(1,:) = FBA_solution.fluxes(excRxns);
    flux(1,:) = FBA_solution.fluxes';
else
    growth_rate = 0;
    exchange_rates(1,:) = zeros(length(excRxns),1);
    flux(1,:) = zeros(size(model.rxns'));
end

% Calculate Biomass and External Metabolite Amounts
biomass(1,:) = biomass_0 + growth_rate.*biomass_0.*dt; % X(t+dt) = X(t) + mu*X(t)*dt
biomass(biomass > max_biomass) = max_biomass;
biomass(biomass < 1E-9) = 0;

% External Metabolite Amounts
mets_conc(1,:) = mets_conc_0' + exchange_rates(1,:).*biomass_0.*dt; % c(t+dt) = c(t) + v*X(t)*dt
mets_conc(mets_conc < 1E-9) = 0; % no negative concentrations

%% Perform Dynamic Flux Balance Analysis - Steps 2 to N

Vuptake_max = Vuptake_max';
for n = 2:N
    % Calculate Flux for Exchange Reactions
    v_uptake = (Vmax.*mets_conc(n-1,:))./(Km + mets_conc(n-1,:)); % v = v_max*([c]/(Km + [c]))
    
    % Check if calculated uptake is greater than max
    idx = find(v_uptake < Vuptake_max);
    v_uptake(idx) = Vuptake_max(idx); % if greater than max, set to max uptake
    v_uptake(v_uptake > -1E-9 & v_uptake < 0) = 0;
    model.lb(excRxns) = v_uptake; % specify new bounds
    
    % FBA
    clear FBA_solution growth_rate
    FBA_solution = FBA(model);
    
    if strcmp(FBA_solution.status,'OPTIMAL') % need optimal solution
        growth_rate = FBA_solution.objectiveValue;
        exchange_rates(n,:) = FBA_solution.fluxes(excRxns);
        flux(n,:) = FBA_solution.fluxes';
    else
        growth_rate = 0;
        exchange_rates(n,:) = zeros(length(excRxns),1);
        flux(n,:) = zeros(size(model.rxns'));
    end
    
    % Calculate Biomass and External Metabolite Amounts
    biomass(n,:) = biomass(n-1,:) + growth_rate.*biomass(n-1,:).*dt; % X(t+dt) = X(t) + mu*X(t)*dt
    biomass(biomass > max_biomass) = max_biomass;
    biomass(biomass < 1E-9) = 0;
    
    % External Metabolite Amounts
    mets_conc(n,:) = mets_conc(n-1,:) + exchange_rates(n,:).*biomass(n-1,:).*dt; % c(t+dt) = c(t) + v*X(t)*dt
    mets_conc(mets_conc < 1E-9) = 0; % no negative concentrations
end
clear n

biomass = [biomass_0; biomass];
mets_conc = [mets_conc_0'; mets_conc];

end