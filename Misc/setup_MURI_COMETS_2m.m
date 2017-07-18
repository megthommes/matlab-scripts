function setup_MURI_COMETS_2m(sparseCon,medium_names,medium_amount,model1,model2,exchMetsFile,forcedSecretion_lb,seed_amt,algSol,model_description,saveLogFolder,dt,N,V,Km,Vmax,saveFilesFolder)
%SETUP_MURI_COMETS_2M Set up models to run in COMETS
%
% setup_MURI_COMETS_2m(sparseCon,medium_names,medium_amount,model1,model2,exchMetsFile,forcedSecretion_lb,seed_amt,algSol,model_description,saveLogFolder)
% setup_MURI_COMETS_2m(sparseCon,medium_names,medium_amount,model1,model2,exchMetsFile,forcedSecretion_lb,seed_amt,algSol,model_description,saveLogFolder,dt,N,V,Km,Vmax,saveFilesFolder)
%
%REQUIRED INPUTS
% medium_names
% medium_amount
% model1,model2
% exchMetsFile
% forcedSecretion_lb
% seed_amt
% algSol
% saveLogFolder: where to save the log files (COMETS Output)
%
%OPTIONAL INPUTS
% model_description
% dt: time step (default = 0.01 hr)
% N: number of steps (default = 500)
% V: volume (default = 100 mL)
% Km (default = 0.01 mM)
% Vmax (default = 10 mmol/gCDW/hr)
% saveFilesFolder (default = current working directory)
%

%% Check Inputs

if (nargin < 9)
    error('myfuns:setup_MURI_COMETS_2m:NotEnoughInputs', ...
        'Not enough input arguments');
end
if ~iscellstr(medium_names)
    error('myfuns:setup_MURI_COMETS_2m:IncorrectInput', ...
        '"medium_names" must be a cell string')
end
if ~isnumeric(medium_amount)
    error('myfuns:setup_MURI_COMETS_2m:IncorrectInput', ...
        '"medium_amount" must be a numeric array')
end
if numel(medium_names) ~= numel(medium_amount)
    error('myfuns:setup_MURI_COMETS_2m:IncorrectInput', ...
        '"medium_names" and "medium_amount" must be the same size')
end
if ~exist('dt','var')
    dt = 0.01; % time step [hr]
end
if ~exist('N','var')
    N = 500; % number of steps
end
if ~exist('V','var')
    V = 100; % volume [mL]
end
if ~exist('Km','var')
    Km = 0.01; % [mM]
end
if ~exist('Vmax','var')
    Vmax = 10; % [mmol/gCDW/hr]
end
if ~exist('saveFilesFolder','var')
    saveFilesFolder = pwd;
end
if isempty(intersect(saveFilesFolder(end),'\'))
    saveFilesFolder = [saveFilesFolder '\'];
end

%% Update Lower Bounds on Exchange Reactions & Make Sure Lower Bound on Biomass is Zero

for tt = 1:numel(sparseCon)
    % Model 1
    model1{tt}.lb(model1{tt}.exch_idx(find(model1{tt}.lb(model1{tt}.exch_idx) == 0))) = -1000;
    model1{tt}.lb(find(model1{tt}.c)) = 0;
    % Model 2
    model2{tt}.lb(model2{tt}.exch_idx(find(model2{tt}.lb(model2{tt}.exch_idx) == 0))) = -1000;
    model2{tt}.lb(find(model2{tt}.c)) = 0;
end

%% Force Secretion & Seed Metabolites

col_letter = num2letter(2*numel(sparseCon));
sparseCon_exchMets = xlsread(exchMetsFile,['A1:' col_letter '1']);
[~,exchType,~] = xlsread(exchMetsFile,['A2:' col_letter '2']);
media_names = arrayfun(@(x) medium_names, 1:numel(sparseCon), 'Uni',false);
media_amount = arrayfun(@(x) medium_amount, 1:numel(sparseCon), 'Uni',false);

% 1->2: Secreted by Model 1
idx_12 = ismember(exchType,'1->2'); idx_12 = find(idx_12 == 1);
for ii = 1:numel(idx_12)
    [~,sparseCon_idx,~] = intersect(sparseCon,sparseCon_exchMets(idx_12(ii)));
    [exchRxns,~] = findExchRxns(model1{sparseCon_idx});
    [~,exchMets,~] = xlsread(exchMetsFile,[num2letter(idx_12(ii)) '3:' num2letter(idx_12(ii)) int2str(sum(exchRxns))]);
    
    % Update Lower Bounds on Forced Secretion
    [~,rxnsIdx,~] = intersect(strrep(strrep(model1{sparseCon_idx}.rxns,'EX_',''),'_e',''),exchMets);
    [~,algSol_rxnsIdx,~] = intersect(algSol.model2{1}.rxns,model1{sparseCon_idx}.rxns(rxnsIdx));
    [~,algSol_sparseConIdx,~] = intersect(algSol.model2{1}.sparse_con,sparseCon(sparseCon_idx));
    if forcedSecretion_lb~=0 && ~isempty(rxnsIdx)
        model1{sparseCon_idx}.lb(rxnsIdx) = forcedSecretion_lb.*algSol.model2{1}.flux(algSol_rxnsIdx,algSol_sparseConIdx);
    end
    
    % Seed Metabolites
    [metsIdx,~] = findMetsFromRxns(model1{sparseCon_idx},rxnsIdx);
    added_media = setdiff(model1{sparseCon_idx}.mets(metsIdx),media_names{sparseCon_idx});
    media_names{sparseCon_idx} = [media_names{sparseCon_idx}; added_media];
    media_amount{sparseCon_idx} = [media_amount{sparseCon_idx}; seed_amt.*ones(size(added_media))];
    
    clear sparseCon_idx exchMets exchRxns rxnsIdx algSol_rxnsIdx added_media metsIdx
end

% 2->1: Secreted by Model 2
idx_21 = ismember(exchType,'2->1'); idx_21 = find(idx_21 == 1);
for ii = 1:numel(idx_21)
    [~,sparseCon_idx,~] = intersect(sparseCon,sparseCon_exchMets(idx_21(ii)));
    [exchRxns,~] = findExchRxns(model2{sparseCon_idx});
    [~,exchMets,~] = xlsread(exchMetsFile,[num2letter(idx_21(ii)) '3:' num2letter(idx_21(ii)) int2str(sum(exchRxns))]);
    
    % Update Lower Bounds on Forced Secretion
    [~,rxnsIdx,~] = intersect(strrep(strrep(model2{sparseCon_idx}.rxns,'EX_',''),'_e',''),exchMets);
    [~,algSol_rxnsIdx,~] = intersect(algSol.model2{2}.rxns,model2{sparseCon_idx}.rxns(rxnsIdx));
    [~,algSol_sparseConIdx,~] = intersect(algSol.model2{2}.sparse_con,sparseCon(sparseCon_idx));
    if forcedSecretion_lb~=0 && ~isempty(rxnsIdx)
        model2{sparseCon_idx}.lb(rxnsIdx) = forcedSecretion_lb.*algSol.model2{2}.flux(algSol_rxnsIdx,algSol_sparseConIdx);
    end
    
    % Seed Metabolites
    [metsIdx,~] = findMetsFromRxns(model2{sparseCon_idx},rxnsIdx);
    added_media = setdiff(model2{sparseCon_idx}.mets(metsIdx),media_names{sparseCon_idx});
    media_names{sparseCon_idx} = [media_names{sparseCon_idx}; added_media];
    media_amount{sparseCon_idx} = [media_amount{sparseCon_idx}; seed_amt.*ones(size(added_media))];
    
    clear sparseCon_idx exchMets exchRxns rxnsIdx algSol_rxnsIdx added_media metsIdx
end

%% Add Model Description

if exist('model_description','var')
    for tt = 1:numel(sparseCon)
        % Model 1
        model1{tt}.description = [model_description ', Model 1, T=' int2str(sparseCon(tt))];
        % Model 2
        model2{tt}.description = [model_description ', Model 2, T=' int2str(sparseCon(tt))];
    end
end

%% Create Layout, Model, & Script Files

for tt = 1:numel(sparseCon)
    model{1} = model1{tt};
    model{2} = model2{tt};    
    % Model 1 Monoculture
    saveFilesPath = [saveFilesFolder 'Model1_Monoculture\T' int2str(sparseCon(tt)) '\'];
    if ~isdir(saveFilesPath); mkdir(saveFilesPath); end
    saveLogPath = [saveLogFolder 'Model1_Monoculture/T' int2str(sparseCon(tt)) '/'];
    createScript(saveFilesPath,[saveLogFolder 'Model1_Monoculture/T' int2str(sparseCon(tt)) '/']);
    world = createLayout(model(1),media_names{tt},media_amount{tt},saveLogPath,dt,N,V,Km,Vmax);
    writeCometsLayout(world,saveFilesPath,'comets_layout.txt',true,false);
    writeCometsModel(model{1},[saveFilesPath getModelName(model{1}) '.txt']);
    % Change Objective Style
    fileID = fopen([saveFilesPath getModelName(model{1}) '.txt'],'a');
    fprintf(fileID,'OBJECTIVE_STYLE\n');
    fprintf(fileID,'\tMAX_OBJECTIVE_MIN_TOTAL\n');
    fprintf(fileID,'//\n');
    fclose(fileID);
    clear world savePath
    % Model 2 Monoculture
    saveFilesPath = [saveFilesFolder 'Model2_Monoculture\T' int2str(sparseCon(tt)) '\'];
    if ~isdir(saveFilesPath); mkdir(saveFilesPath); end
    saveLogPath = [saveLogFolder 'Model2_Monoculture/T' int2str(sparseCon(tt)) '/'];
    createScript(saveFilesPath,[saveLogFolder 'Model2_Monoculture/T' int2str(sparseCon(tt)) '/']);
    world = createLayout(model(2),media_names{tt},media_amount{tt},saveLogPath,dt,N,V,Km,Vmax);
    writeCometsLayout(world,saveFilesPath,'comets_layout.txt',true,false);
    writeCometsModel(model{2},[saveFilesPath getModelName(model{2}) '.txt']);
    % Change Objective Style
    fileID = fopen([saveFilesPath getModelName(model{2}) '.txt'],'a');
    fprintf(fileID,'OBJECTIVE_STYLE\n');
    fprintf(fileID,'\tMAX_OBJECTIVE_MIN_TOTAL\n');
    fprintf(fileID,'//\n');
    fclose(fileID);
    clear world savePath
    % Coculture
    saveFilesPath = [saveFilesFolder 'Coculture\T' int2str(sparseCon(tt)) '\'];
    if ~isdir(saveFilesPath); mkdir(saveFilesPath); end
    saveLogPath = [saveLogFolder 'Coculture/T' int2str(sparseCon(tt)) '/'];
    createScript(saveFilesPath,[saveLogFolder 'Coculture/T' int2str(sparseCon(tt)) '/']);
    world = createLayout(model,media_names{tt},media_amount{tt},saveLogPath,dt,N,V,Km,Vmax);
    writeCometsLayout(world,saveFilesPath,'comets_layout.txt',true,false);
    for num_model = 1:2
        writeCometsModel(model{num_model},[saveFilesPath getModelName(model{num_model}) '.txt']);
        % Change Objective Style
        fileID = fopen([saveFilesPath getModelName(model{num_model}) '.txt'],'a');
        fprintf(fileID,'OBJECTIVE_STYLE\n');
        fprintf(fileID,'\tMAX_OBJECTIVE_MIN_TOTAL\n');
        fprintf(fileID,'//\n');
        fclose(fileID);
    end
    clear world savePath
end

end

% COMETS Layout
function [world] = createLayout(model,medium_names,medium_amount,saveLogFolder,dt,N,V,Km,Vmax)
% world = createLayout(model,medium_names,medium_amount,saveLogFolder)
% world = createLayout(model,medium_names,medium_amount,saveLogFolder,dt,N,V)
%
%REQUIRED INPUTS
% model
% medium_names
% medium_amount
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
world.params.timeStep = dt; % time step [hr]
world.params.maxCycles = N; % number of steps

% Set Initial Medium
for met_num = 1:numel(medium_names)
    world = world.setInitialMedia(medium_names(met_num),medium_amount(met_num));
end

% Volume
world.params.spaceWidth = V^(1/3); % [cm]
% volume is spaceWidth^3 (1 mL = 1 cm^3)

% Make Sure Death Rate is Zero
world.params.deathRate = 0;

% Make Sure Cell Overlap is Allowed
world.params.allowCellOverlap = true;

% Set Km
world.params.defaultKm = Km*1e-3; % [mmol/cm^3]
% [mmol/cm^3] -> mM = mmol/L * 1 L/1e3 mL * 1 mL/1 cm^3

% Set Vmax
world.params.defaultVmax = Vmax; % [mmol/gCDW/hr]

% Set Maximum Biomass
world.params.maxSpaceBiomass = 1*V*1e-3; % [gCDW] -> 1 g/L * 1 L *1 L/1e-3 mL
% http://bionumbers.hms.harvard.edu/bionumber.aspx?&id=104943&ver=15&trm=e.%20coli%20inoculum%20concentration

end

% COMETS Script File
function createScript(directory,layoutPath)
% Modified from createScriptFile
    scriptfile = fullfile(directory, 'comets_script.txt'); % path to script file
    lfile = fullfile(layoutPath,'comets_layout.txt'); % path to layout file
    lfile = strrep(lfile,'\','/');
    file = fopen(scriptfile,'w'); % open script file
    fprintf(file,'load_layout %s\n',lfile);
    fclose(file); % close script file
end







