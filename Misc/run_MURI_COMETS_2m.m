function [time,biomass,flux,metAmount,metNames,cometsOutput] = run_MURI_COMETS_2m(model,medium_names,medium_amount,saveLogFolder,dt,N)
%RUN_MURI_COMETS_2M Run models in COMETS
%
% [time,biomass,flux,metAmount,metNames,cometsOutput] = run_MURI_COMETS_2m(model,medium_names,medium_amount,saveLogFolder)
% [time,biomass,flux,metAmount,metNames,cometsOutput] = run_MURI_COMETS_2m(model,medium_names,medium_amount,saveLogFolder,dt,N)
%
%REQUIRED INPUTS
% model
% medium_names
% medium_amount
% saveLogFolder
%
%OPTIONAL INPUT
% dt: time step (default = 0.01)
% N: number of steps (default = 200)
%
%OUTPUTS
% time
% biomass
% flux
% metAmount
% metNames
% cometsOutput

%% Check Inputs

if (nargin < 4)
    error('myfuns:run_MURI_COMETS_2m:NotEnoughInputs', ...
        'Not enough input arguments');
end
if ~isdir(saveLogFolder)
    mkdir(saveLogFolder); % create folder if it does not already exist
end

cwd = pwd;

%% COMETS Layout

% Create COMETS Layout
world = CometsLayout();

% Add Models
for num_model = 1:numel(model) % add model(s)
    world = world.addModel(model{num_model});
end

% Write Logs
world.params.writeBiomassLog = true; % write biomass log
world.params.biomassLogName = [saveLogFolder world.params.biomassLogName(3:end)]; % save in specific folder
world.params.writeMediaLog = true; % write media log
world.params.mediaLogName = [saveLogFolder world.params.mediaLogName(3:end)]; % save in specific folder
world.params.writeFluxLog = true; % write flux log
world.params.fluxLogName = [saveLogFolder world.params.fluxLogName(3:end)]; % save in specific folder

% Initialize Population to Default
world = setInitialPop(world,'1x1'); % initialize population to default (1E-5)

% Time
if exist('dt','var')
    world.params.timeStep = dt; % time step [hr]
end
if exist('N','var')
    world.params.maxCycles = N; % number of steps
end

% Set Initial Medium
for met_num = 1:numel(medium_names)
    world = world.setInitialMedia(medium_names(met_num),medium_amount(met_num));
end

% Volume
world.params.spaceWidth = 10; % [cm]
% volume is spaceWidth^3 (1 L = 1E3 cm^3 = (10 cm)^3)

% Make Sure Death Rate is Zero
world.params.deathRate = 0;

% Make Sure Cell Overlap is Allowed
world.params.allowCellOverlap = true;

%% Run Simulation

cd(saveLogFolder)
% cometsOutput = runComets(world);
createCometsFiles(world,saveLogFolder);%create layout,model,script & param files
[~,comets_path] = system('echo %COMETS_HOME%'); % determine COMETS path
comets_path = strtrim(comets_path); % remove leading and trailing white space
javaclasspath(comets_path); % add COMETS classpath to MATLAB
fid = fopen([saveLogFolder '\comets_w64_scr.bat'],'w');
fprintf(fid,'%s%s%s','java -Xmx2048m -classpath %COMETS_HOME%\bin\COMETS_2.0.3.jar;%COMETS_HOME%\lib\x64\glpk-java.jar;%COMETS_HOME%\lib\jogamp-all-platforms\jar\jogl-all.jar;%COMETS_HOME%\lib\jogamp-all-platforms\jar\gluegen.jar;%COMETS_HOME%\lib\jogamp-all-platforms\jar\gluegen-rt.jar;%COMETS_HOME%\lib\jogamp-all-platforms\jar\gluegen-rt-natives-windows-amd64.jar;%COMETS_HOME%\lib\jogamp-all-platforms\jar\jogl-all-natives-windows-amd64.jar;%GUROBI_HOME%\lib\gurobi.jar    -Djava.library.path=%COMETS_HOME%\lib\x64;%GUROBI_HOME%\lib\gurobi56.lib;%GUROBI_HOME%\bin  edu.bu.segrelab.comets.Comets -loader edu.bu.segrelab.comets.fba.FBACometsLoader -script ', saveLogFolder, 'comets_script.txt');
fclose(fid);
[status,cometsOutput] = system('comets_w64_scr.bat');
delete([saveLogFolder 'comets_w64_scr.bat']); % remove script bat file from folder
cd(cwd)

% Save COMETS Output
fid = fopen([saveLogFolder '\comets_output.txt'],'w');
fprintf(fid,'%s',cometsOutput);
fclose(fid);

%% Parse Logs

% if status == 0
% Biomass
biomass_out = parseBiomassCometsOutput(world.params.biomassLogName);
biomass = cell(numel(model),1);
for num_model = 1:numel(model)
    biomass{num_model} = biomass_out.biomass(biomass_out.model == num_model-1);
end

% Time
N = numel(biomass{1})-1;
time = 0:dt:N*dt;

% Medium
media_out = parseMediaLog(world.params.mediaLogName);
[metNames,~,met_idx] = unique(cellfun(@char,media_out.metname, 'Uni',false));
if numel(media_out.amt)/numel(metNames) ~= N+1 % Check Size of Medium
    idx = find(media_out.t > N+1);
    media_out(idx,:) = []; % remove extra time steps
    met_idx(idx) = [];
end
metAmount = zeros(N+1,numel(metNames));
for num_met = 1:numel(metNames)
    metAmount(:,num_met) = media_out.amt(met_idx == num_met);
end

% Flux
flux_out = parseFluxLog(world.params.fluxLogName, world);
flux = cell(numel(model),1);
for num_model = 1:numel(model)
    flux{num_model} = cell2mat(arrayfun(@(t) flux_out.flux(flux_out.model == num_model & flux_out.t == t),1:N, 'Uni',false))';
end

fclose('all');

end