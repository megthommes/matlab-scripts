function [output] = growth_vs_uptake(model,carbonRxn,nitrogenRxn,folder)
%growth_vs_uptake Generates growth vs. uptake rate curves
%   Generates phenotype phase planes
%
%
%
%   Meghan Thommes
%   11/24/14

initCobraToolbox

%% Check Inputs

if nargin < 2
    error('Not enough parameters');
elseif nargin == 3
    if isempty(nitrogenRxn) || ~ischar(nitrogenRxn)
        N_flag = 0;
    else
        N_flag = 1;
    end
elseif nargin == 4
    if isempty(nitrogenRxn) || ~ischar(nitrogenRxn)
        N_flag = 0;
    else
        N_flag = 1;
    end
    
    if isempty(folder) || ~ischar(folder)
        folder = pwd;
    elseif exist(folder,'file') ~= 7
        mkdir(folder);
    end
end

%% Vary Carbon

% Carbon Exchange Reaction
carbonRxn_idx = strcmp(model.rxns,carbonRxn);

% Carbon Name
temp = model.S(:,carbonRxn_idx);
[carbon_ind,~] = find(temp ~= 0);
carbonName = model.metNames(carbon_ind);
carbonName = strrep(carbonName,'_',' ');
carbonName = carbonName{:};
clear temp carbon_ind

lb = -1.*[1 10 25 50 75 100:100:1000];
growthRate_C = zeros(size(lb));
carbonRate_C = zeros(size(lb));

for C = 1:length(lb)
    model_temp = changeRxnBounds(model,carbonRxn,lb(C),'l');
    solution = optimizeCbModel(model_temp);
    
    % Growth Rate [mmol/(gDW*hr)]
    growthRate_C(C) = solution.f;
    
    if solution.stat ~= 0 % not infeasible
        % Carbon Update Rate [mmol/(gDW*hr)]
        carbonRate_C(C) = model_temp.lb(carbonRxn_idx); % uptake rate = - lower bound
        carbonRxnRate_C(C) = solution.x(carbonRxn_idx); % flux through reaction
        
        % Flux Solution
        fluxes_C(:,C) = solution.x;
    else
        % Carbon Update Rate [mmol/(gDW*hr)]
        carbonRate_C(C) = NaN; % uptake rate = - lower bound
        carbonRxnRate_C(C) = NaN; % flux through reaction
        
        % Flux Solution
        fluxes_C(:,C) = NaN;
    end
    clear model_temp solution
end

%% Vary Nitrogen

if N_flag == 1
    % Nitrogen Exchange Reaction
    nitrogenRxn_idx = strcmp(model.rxns,nitrogenRxn);
    
    % Nitrogen Name
    temp = model.S(:,nitrogenRxn_idx);
    [nitrogen_ind,~] = find(temp ~= 0);
    nitrogenName = model.metNames(nitrogen_ind);
    nitrogenName = strrep(nitrogenName,'_',' ');
    nitrogenName = nitrogenName{:};
    clear temp nitrogen_ind
    
    lb = -1.*[1 10 25 50 75 100:100:1000];
    growthRate_N = zeros(length(lb),length(lb));
    carbonRate_N = zeros(length(lb),length(lb));
    nitrogenRate_N = zeros(length(lb),length(lb));
    
    for C = 1:length(lb) % carbon
        for N = 1:length(lb) % nitrogen
            model_temp = changeRxnBounds(model,nitrogenRxn,lb(N),'l');
            model_temp = changeRxnBounds(model_temp,carbonRxn,lb(C),'l');
            solution = optimizeCbModel(model_temp);
            
            % Growth Rate [mmol/(gDW*hr)]
            growthRate_N(C,N) = solution.f;
            
            if solution.stat ~= 0 % not infeasible
                % Carbon Update Rate [mmol/(gDW*hr)]
                carbonRate_N(C,N) = model_temp.lb(carbonRxn_idx); % uptake rate = - lower bound
                nitrogenRate_N(C,N) = model_temp.lb(nitrogenRxn_idx); % uptake rate = - lower bound
                
                % Flux Solution
                fluxes_N(:,C,N) = solution.f;
            else
                % Carbon Update Rate [mmol/(gDW*hr)]
                carbonRate_N(C,N) = NaN; % uptake rate = - lower bound
                nitrogenRate_N(C,N) = NaN; % uptake rate = - lower bound
                
                % Flux Solution
                fluxes_N(:,C,N) = NaN;
            end          
            
            clear model_temp solution
        end
    end
end

%% Plot

min_lb = carbonRate_C(find(growthRate_C == max(growthRate_C),1)); % min lower bound to get max growth rate

figC = figure(1); % new figure
plot(-carbonRate_C,growthRate_C,'rx--','MarkerSize',8,'LineWidth',2);
xlabel([carbonName ' Uptake Rate [mmol/(gDW \cdot hr)]'],'fontsize',12)
ylabel('Growth Rate [mmol/(gDW \cdot hr)]','fontsize',12)
title(strrep(model.description,'_',' '),'fontsize',16)
text('Position',[-min_lb max(growthRate_C)-2], ...
    'String',['Max Growth Rate = ' num2str(max(growthRate_C))], ...
    'BackgroundColor','w')

if N_flag == 1
    figN = figure(2); % new figure
    plot(-carbonRate_N,growthRate_N,'x--','MarkerSize',8,'LineWidth',2);
    h = legend(cellstr(num2str((nitrogenRate_N(7,:)'))),'Location','EastOutside');
    set(get(h,'title'),'string',{nitrogenName, 'Uptake Rate'}); % legend title
    set(h,'EdgeColor','w'); % remove box
    xlabel([strrep(carbonName,'_',' ') ' Uptake Rate [mmol/(gDW \cdot hr)]'],'fontsize',12)
    ylabel('Growth Rate [mmol/(gDW \cdot hr)]','fontsize',12)
    title(strrep(model.description,'_',' '),'fontsize',16)
end

%% Save Variables, Figures, and Data

output.lb = lb;
output.carbonRxn = carbonRxn;
output.carbonName = carbonName;
output.growthRate_C = growthRate_C;
output.carbonRate_C = carbonRate_C;
output.carbonRxnRate_C = carbonRxnRate_C;
output.fluxes_C = fluxes_C;
output.growth_figC = figC;
saveas(figC,[folder '\growth_varyC'],'png')

if N_flag == 1
    output.nitrogenRxn = nitrogenRxn;
    output.nitrogenName = nitrogenName;
    output.growthRate_N = growthRate_N;
    output.carbonRate_N = carbonRate_N;
    output.nitrogenRate_N = nitrogenRate_N;
    output.fluxes_N = fluxes_N;
    output.growth_figN = figN;
    saveas(figN,[folder '\growth_varyN'],'png')
end

save([folder '\growthCurves.mat'],'-struct','output');
modelName = model.description;
save([folder '\growthCurves.mat'],'modelName','-append')

end






