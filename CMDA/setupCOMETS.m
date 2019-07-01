function setupCOMETS(fba_models,medium,sec_lb,seed_amt,params,saveLogFolder,saveFilesFolder)
%setupCOMETS Create COMETS model and layouts to run
%
%setupCOMETS(fba_models,medium)
%setupCOMETS(fba_models,medium,sec_lb,seed_amt,params,saveLogFolder,saveFilesFolder)
%
%REQUIRED INPUTS
% fba_models: structure with fields:
%   S: Stoichiometric matrix
%   b: Right hand side = dx/dt
%   c: Objective coefficients
%   lb: Lower bounds
%   ub: Upper bounds
%   vmax: vector, size of exch_rxns
%       **Optional Field (default = 10 mmol/gCDW/hr)
%   km: vector, size of exch_rxns
%       **Optional Field (default = 0.01 mM (1E-5 mmol/cm^3))
% medium: Structure with fields:
%   names: names of medium metabolites
%   amount: amount of medium metabolites (mmol)
%   static: logical array indicating where medium is static (1: static)
%       **Optional Field (default: all 0)
%
%OPTIONAL INPUTS
% sec_lb: forced secretion, fraction of algorithm flux value (default = 0)
% seed_amt: exchanged metabolites seed amount (default = 0 mmol)
% params: Structure with fields:
%   dt: time step (default = 0.01 hr)
%   N: number of steps (default = 500)
%   V: volume (default = 100 mL)
% saveLogFolder: where to save the COMETS output log files (default = current working directory)
% saveFilesFolder: directory to save files (default = current working directory)
%
% Meghan Thommes 02/05/2018

%% Check Inputs

% Required
if (nargin < 2)
    error('myfuns:setupCOMETS:NotEnoughInputs', ...
        'Not enough inputs: need fba_models and medium');
else
    % fba_models
    if ~iscell(fba_models)
        error('myfuns:setupCOMETS:IncorrectType', ...
            '"fba_models" needs to be a cell structure');
    else
        numModels = numel(fba_models);
        for model_num = 1:numModels
            if ~isfield(fba_models{model_num},'vmax')
                fba_models{model_num}.vmax = 10.*ones(size(fba_models{model_num}.exch_idx));
            end
            if ~isfield(fba_models{model_num},'km')
                fba_models{model_num}.km = 0.01.*1e-3.*ones(size(fba_models{model_num}.exch_idx));
            end
        end
    end
    % medium
    if ~isstruct(medium)
        error('myfuns:setupCOMETS:IncorrectType', ...
            '"medium" needs to be a structure');
    elseif ~isfield(medium,'names') || ~isfield(medium,'amount')
        error('myfuns:setupCOMETS:IncorrectType', ...
            '"medium" needs "names" and "amount" fields');
    end
end

% Optional
if ~isfield(medium,'static')
    medium.static = logical(zeros(size(medium.names)));
end
if ~exist('sec_lb','var')
    sec_lb = 0; % forced secretion
end
if ~exist('seed_amt','var')
    seed_amt = 0; % exchanged metabolites seed amount
end
if ~isfield(params,'dt')
    params.dt = 0.01; % time step [hr]
end
if ~isfield(params,'N')
    params.N = 500; % number of steps
end
if ~isfield(params,'V')
    params.V = 100; % volume [mL]
end
if ~exist('saveFilesFolder','var')
    saveFilesFolder = pwd;
end

% "Default" Vmax & Km
Vmax = []; Km = [];
for model_num = 1:numModels
    Vmax = [Vmax; fba_models{model_num}.vmax];
    Km = [Km; fba_models{model_num}.km];
end
% Vmax
uni_Vmax = unique(Vmax);
N = histcounts(Vmax,numel(uni_Vmax));
params.Vmax = uni_Vmax(find(N == max(N)));
% Km
uni_Km = unique(Km);
N = histcounts(Km,numel(uni_Km));
params.Km = uni_Km(find(N == max(N)));
clear Vmax Km uni_* N

%% Update Lower Bounds on Exchange Reactions & Make Sure Lower Bound on Biomass is Zero

for model_num = 1:numModels
    % Update Lower Bounds on Exchange Reactions
    idx = find(fba_models{model_num}.lb(fba_models{model_num}.exch_idx) == 0);
    fba_models{model_num}.lb(fba_models{model_num}.exch_idx(idx)) = min(fba_models{model_num}.lb(fba_models{model_num}.exch_idx));
    
    % Make Sure Lower Bound on Biomass is Zero
    fba_models{model_num}.lb(fba_models{model_num}.bio_idx) = 0;
end

%% Force Secretion & Seed Metabolites

if sec_lb > 0 % if have forced secretion
    % Extracellular Metabolites
    [exchMets_idx,~] = identifyExchMets(fba_models{1},fba_models{1}.exch_idx);
    
    % Exchanged Metabolites
    exchOrder = combnk(1:numModels,2);
    for ii = 1:size(exchOrder,1)
        m1 = exchOrder(ii,1);
        m2 = exchOrder(ii,2);
        
        % m1 -> m2
        fieldName = ['exchMets_' int2str(m1) int2str(m2)];
        if isfield(fba_models{m1},fieldName)
            [~,~,exchMetsIdx_m1m2] = intersect(fba_models{m1}.(fieldName),fba_models{m1}.metNames(exchMets_idx),'stable');
            fba_models{m1}.lb(fba_models{m1}.exch_idx(exchMetsIdx_m1m2)) = sec_lb.*fba_models{m1}.flux(fba_models{m1}.exch_idx(exchMetsIdx_m1m2));
            % update medium
            clear idx
            [~,idx] = setdiff(fba_models{m1}.mets(exchMetsIdx_m1m2), medium.names);
            medium.names = [medium.names; fba_models{m1}.mets(exchMetsIdx_m1m2(idx))];
            medium.amount = [medium.amount; seed_amt.*ones(numel(idx),1)];
        end
        
        % m2 -> m1
        fieldName = ['exchMets_' int2str(m2) int2str(m1)];
        if isfield(fba_models{m2},fieldName)
            [~,~,exchMetsIdx_m2m1] = intersect(fba_models{m2}.(fieldName),fba_models{m2}.metNames(exchMets_idx),'stable');
            fba_models{m2}.lb(fba_models{m2}.exch_idx(exchMetsIdx_m2m1)) = sec_lb.*fba_models{m2}.flux(fba_models{m2}.exch_idx(exchMetsIdx_m2m1));
            % update medium
            clear idx
            [~,idx] = setdiff(fba_models{m2}.mets(exchMetsIdx_m2m1), medium.names);
            medium.names = [medium.names; fba_models{m2}.mets(exchMetsIdx_m2m1(idx))];
            medium.amount = [medium.amount; seed_amt.*ones(numel(idx),1)];
        end
    end
else % if don't have forced secretion
    % Extracellular Metabolites
    [exchMets_idx,~] = identifyExchMets(fba_models{1},fba_models{1}.exch_idx);
    
    % Exchanged Metabolites
    exchOrder = combnk(1:numModels,2);
    for ii = 1:size(exchOrder,1)
        m1 = exchOrder(ii,1);
        m2 = exchOrder(ii,2);
        
        % m1 -> m2
        fieldName = ['exchMets_' int2str(m1) int2str(m2)];
        if isfield(fba_models{m1},fieldName)
            [~,~,exchMetsIdx_m1m2] = intersect(fba_models{m1}.(fieldName),fba_models{m1}.metNames(exchMets_idx),'stable');
            % update medium
            clear idx
            [~,idx] = setdiff(fba_models{m1}.mets(exchMetsIdx_m1m2), medium.names);
            medium.names = [medium.names; fba_models{m1}.mets(exchMetsIdx_m1m2(idx))];
            medium.amount = [medium.amount; seed_amt.*ones(numel(idx),1)];
        end
        
        % m2 -> m1
        fieldName = ['exchMets_' int2str(m2) int2str(m1)];
        if isfield(fba_models{m2},fieldName)
            [~,~,exchMetsIdx_m2m1] = intersect(fba_models{m2}.(fieldName),fba_models{m2}.metNames(exchMets_idx),'stable');
            % update medium
            clear idx
            [~,idx] = setdiff(fba_models{m2}.mets(exchMetsIdx_m2m1), medium.names);
            medium.names = [medium.names; fba_models{m2}.mets(exchMetsIdx_m2m1(idx))];
            medium.amount = [medium.amount; seed_amt.*ones(numel(idx),1)];
        end
    end
end

%% Create Layout, Model, & Script Files

% Override Default Kinetic Parameters
cometsParams = CometsParams();
cometsParams.defaultReactionUpper =  1000;
cometsParams.defaultReactionLower = -1000;
cometsParams.defaultKm = params.Km;
cometsParams.defaultVmax = params.Vmax;
cometsParams.objectiveStyle = 'MAX_OBJECTIVE_MIN_TOTAL';

% Monoculture
for model_num = 1:numModels
    % File Paths
    simFolder = ['Model' int2str(model_num) '_Monoculture\'];
    saveFilesPath = [saveFilesFolder simFolder];
    if ~isdir(saveFilesPath); mkdir(saveFilesPath); end
    saveLogPath = [saveLogFolder simFolder]; saveLogPath = strrep(saveLogPath,'\','/');
    % Script File
    createScript(saveFilesPath,[saveLogFolder simFolder]);
    % Layout File
    model{1} = fba_models{model_num};
    world = makeLayout(model,medium,params,saveLogPath); clear model
    writeCometsLayout(world,saveFilesPath,'comets_layout.txt',true,false);
    % Model File
    writeCometsModel(fba_models{model_num},[saveFilesPath getModelName(fba_models{model_num}) '.txt'],cometsParams);
    clear world
end

% Coculture
if numModels > 1
    % File Paths
    simFolder = 'Coculture\';
    saveFilesPath = [saveFilesFolder simFolder];
    if ~isdir(saveFilesPath); mkdir(saveFilesPath); end
    saveLogPath = [saveLogFolder simFolder]; saveLogPath = strrep(saveLogPath,'\','/');
    % Script File
    createScript(saveFilesPath,[saveLogFolder simFolder]);
    % Layout File
    model = cell(numModels,1);
    for model_num = 1:numModels
        model{model_num} = fba_models{model_num};
    end
    world = makeLayout(model,medium,params,saveLogPath); clear model
    writeCometsLayout(world,saveFilesPath,'comets_layout.txt',true,false);
    % Model File
    for model_num = 1:numModels
        writeCometsModel(fba_models{model_num},[saveFilesPath getModelName(fba_models{model_num}) '.txt'],cometsParams);
    end
    clear world
end

end

%% COMETS Script File
function createScript(directory,layoutPath)
% Modified from createScriptFile
scriptfile = fullfile(directory, 'comets_script.txt'); % path to script file
lfile = fullfile(layoutPath,'comets_layout.txt'); % path to layout file
lfile = strrep(lfile,'\','/');
file = fopen(scriptfile,'w'); % open script file
fprintf(file,'load_layout %s\n',lfile);
fclose(file); % close script file
end

%% COMETS Layout
function [world] = makeLayout(model,medium,params,saveLogFolder)
% world = makeLayout(model,medium,params,saveLogFolder)
%
%REQUIRED INPUTS
% model: structure with fields:
%   S: Stoichiometric matrix
%   b: Right hand side = dx/dt
%   c: Objective coefficients
%   lb: Lower bounds
%   ub: Upper bounds
% medium: Structure with fields:
%   names: names of medium metabolites
%   amount: amount of medium metabolites (mmol)
%   static: logical array indicating where medium is static (1: static)
%       **Optional Field (default: all 0)
% params: Structure with fields:
%   dt: time step
%   N: number of steps
%   V: volume
%   Vmax: maximum uptake rate of exchange reactions
%   Km: 
% saveLogFolder

% Create COMETS Layout
world = CometsLayout();

% Add Models
for num_model = 1:numel(model)
    world = world.addModel(model{num_model});
end

% Write Logs
world.params.writeBiomassLog = true; % write biomass log
world.params.biomassLogName = [saveLogFolder world.params.biomassLogName(3:end)]; % save in specific folder
world.params.writeMediaLog = true; % write media log
world.params.mediaLogName = [saveLogFolder world.params.mediaLogName(3:end)]; % save in specific folder
world.params.writeFluxLog = true; % write flux log
world.params.fluxLogName = [saveLogFolder world.params.fluxLogName(3:end)]; % save in specific folder
world.params.writeMatFile = true; % write everything
world.params.MatFileName = [saveLogFolder world.params.MatFileName(3:end)]; % save in specific folder

% Initialize Population to Default
world = setInitialPop(world,'1x1'); % initialize population to default
world.initial_pop = 1E-7.*ones(numel(model),1);

% Time
world.params.timeStep = params.dt; % time step [hr]
world.params.maxCycles = params.N; % number of steps

% Define Medium
for met_num = 1:numel(medium.names)
    world = world.setInitialMedia(medium.names(met_num),medium.amount(met_num));
end

% Set Static Reactions
[~,~,idx] = intersect(medium.names(logical(medium.static)),world.mets,'stable');
world.global_static_media(idx,1) = 1;
world.global_static_media(idx,2) = medium.amount(logical(medium.static));

% Volume
world.params.spaceWidth = params.V^(1/3); % [cm]
% volume is spaceWidth^3 (1 mL = 1 cm^3)

% Make Sure Death Rate is Zero
world.params.deathRate = 0;

% Make Sure Cell Overlap is Allowed
world.params.allowCellOverlap = true;

% Set Km
world.params.defaultKm = params.Km; % [mmol/cm^3]

% Set Vmax
world.params.defaultVmax = params.Vmax; % [mmol/gCDW/hr]

% Set Maximum Biomass
world.params.maxSpaceBiomass = 1*params.V*1e-3; % [gCDW] -> 1 g/L * 1 L *1 L/1e-3 mL
% http://bionumbers.hms.harvard.edu/bionumber.aspx?&id=104943&ver=15&trm=e.%20coli%20inoculum%20concentration

end
    
    