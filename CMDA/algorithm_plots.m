function fig_handles = algorithm_plots(model_flux1,model_flux2,model_title,plotFolder)

%INPUTS
% model_flux*: cell structure that contains the following fields: (fields
%  that are only present when there are 2 models are marked with *)
%   intl_con: intracellular sparsity constraints, vector
%   trspt_con: transport sparsity constraints, vector
%   mets: metabolite names, vector
%   rxns: reaction names, vector
%   biomass: biomass flux, vector
%   flux: reaction fluxes, matrix [rxns x sparse_con]
%   int: reaction binary variables t, matrix [rxns x sparse_con]
%   exch_idx: index of exchange ("uptake") reactions
%   trspt_idx: index of transport ("exchange") reactions
%   intl_idx: index of intracellular ("internal") reactions
%   mediumMets: names of metabolites in the medium
%   mediumRxns: names of reactions in the medium
%   warningFlag: 0 if no warning, 1 if warning
%       BinaryVariableFlux: Constraint violation, model has non-zero flux when t_i=0
%       CalcExchFlux: Incorrect calculation of exchange flux
%   exchMets: metabolite names of calc/totalExch_flux
%   totalExch_flux: total exchange flux from CMDA
%   totalCalcExch_flux: calculated total exchange flux
%   fluxType: Number indicates what type of flux occurs [time x metabolites]
%       0: not used by model 1 or model 2 (->)
%       1: produced by model 1 and model 2 (1,2->)
%       2: produced by model 1 and not used by model 2 (1->)
%       3: not used by model 1 and produced by model 2 (2->)
%       4: consumed by model 1 and model 2 (->1,2)
%       5: consumed by model 1 and not used by model 2 (->1)
%       6: not used by model 1 and consumed by model 2 (->2)
%       7: consumed by model 1 and produced by model 2 (2->1)
%       8: produced by model 1 and consumed by model 2 (1->2)
%   exchangedMetabolites
%       summary: rows are the metabolite names, each column shows the
%           number of times a metabolite was exchanged (1->2, 2->1, 1->2 & 2->1)
%       model1_to_model2/model2_to_model1: metabolites exchanged at each
%           intracellular sparsity constraint
% model_title: name of model
% plotFolder: cell array path to save figures. First index is for
%  model_flux1, second index is for model_flux2
%   If not included, figures are not saved.
%
% If you are not plotting both model_flux or if you are not saving the
%  figure, then enter "''" for the model_flux you are not
%  plotting/plotFolder.
%   Example: algorithm_plots(model_flux1,model_flux2,plotFolder,model_title)
%
%OUTPUT
% fig_handles
%
% Meghan Thommes 08/15/2017

%% Check Inputs

if ~isempty(plotFolder)
    saveFlag = 1;
    if ~isempty(model_flux1) && ~isempty(model_flux2)
        plotFolder1 = plotFolder{1};
        plotFolder2 = plotFolder{2};
        intlCon_idx1 = find(model_flux1{1}.biomass~=0);
        intlCon_idx2 = find(model_flux2{1}.biomass~=0);
        plotFlag = 3;
    elseif isempty(model_flux1)
        plotFolder2 = plotFolder{2};
        intlCon_idx2 = find(model_flux2{1}.biomass~=0);
        plotFlag = 2;
    elseif isempty(model_flux2)
        plotFolder1 = plotFolder{1};
        intlCon_idx1 = find(model_flux1{1}.biomass~=0);
        plotFlag = 1;
    end
else
    saveFlag = 0;
    if ~isempty(model_flux1) && ~isempty(model_flux2)
        intlCon_idx1 = find(model_flux1{1}.biomass~=0);
        intlCon_idx2 = find(model_flux2{1}.biomass~=0);
        plotFlag = 3;
    elseif isempty(model_flux1)
        intlCon_idx2 = find(model_flux2{1}.biomass~=0);
        plotFlag = 2;
    elseif isempty(model_flux2)
        intlCon_idx1 = find(model_flux1{1}.biomass~=0);
        plotFlag = 1;
    end
end

%% Variables

% Plots (because MATLAB isn't using the default values)
legendSize = 12;
axesLabelSize = 14;
xyLabelSize = 18;
titleSize = 22;
lineWidth = 2;
markerSize = 75;
tickLength = 0.005;

% Biomass Colors
color_1m = [9 110 56]./256; % dark green
color_2m1 = [200 30 41]./256; % dark red
color_2m2 = [76 80 256]./256; % dark blue/purple

%% Biomass Flux

if plotFlag == 3 % 1 Model & 2 Models
    % All Intracellular Sparsity Constraints
    fig_handles(1) = figure(1);
    plot(model_flux1{1}.intl_con,model_flux1{1}.biomass,'s-', 'Color',color_1m, 'LineWidth',lineWidth); hold on
    plot(model_flux2{1}.intl_con,model_flux2{1}.biomass,'^--', 'Color',color_2m1, 'LineWidth',lineWidth);
    plot(model_flux2{1}.intl_con,model_flux2{1}.biomass,'v:', 'Color',color_2m2, 'LineWidth',lineWidth); hold off
    xlim([min([model_flux1{1}.intl_con, model_flux2{1}.intl_con]) max([model_flux1{1}.intl_con, model_flux2{1}.intl_con])])
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux'], 'FontSize',titleSize)
    legend({'1 Model','2 Models, Model #1','2 Models, Model #2'}, ...
        'FontSize',legendSize, 'Location','SouthEast', 'EdgeColor','w');
    
    % Only Intracellular Sparsity Constraints With Growth and No Plateau
    model1_max = model_flux1{1}.intl_con(intlCon_idx1(find(abs(diff(model_flux1{1}.biomass(intlCon_idx1))) <= 1E-9,1))); if isempty(model1_max); model1_max = model_flux1{1}.intl_con(end); end
    model2_max = model_flux2{1}.intl_con(intlCon_idx2(find(abs(diff(model_flux2{1}.biomass(intlCon_idx2))) <= 1E-9,1))); if isempty(model2_max); model2_max = model_flux2{1}.intl_con(end); end
    model1_min = model_flux1{1}.intl_con(find(model_flux1{1}.biomass == 0, 1, 'last')); if isempty(model1_min); model1_min = model_flux1{1}.intl_con(1); end
    model2_min = model_flux2{1}.intl_con(find(model_flux2{1}.biomass == 0, 1, 'last')); if isempty(model2_min); model2_min = model_flux2{1}.intl_con(1); end
    fig_handles(2) = figure(2);
    plot(model_flux1{1}.intl_con,model_flux1{1}.biomass,'s-', 'Color',color_1m, 'LineWidth',lineWidth); hold on
    plot(model_flux2{1}.intl_con,model_flux2{1}.biomass,'^--', 'Color',color_2m1, 'LineWidth',lineWidth);
    plot(model_flux2{1}.intl_con,model_flux2{1}.biomass,'v:', 'Color',color_2m2, 'LineWidth',lineWidth); hold off
    xlim([nanmin(model1_min,model2_min)-1 nanmax(model1_max,model2_max)+1])
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux'], 'FontSize',titleSize)
    legend({'1 Model','2 Models, Model #1','2 Models, Model #2'}, ...
        'FontSize',legendSize, 'Location','SouthEast', 'EdgeColor','w');
    
elseif plotFlag == 2 % 2 Models
    % All Intracellular Sparsity Constraints
    fig_handles(1) = figure(1);
    plot(model_flux2{1}.intl_con,model_flux2{1}.biomass,'^--', 'Color',color_2m1, 'LineWidth',lineWidth); hold on
    plot(model_flux2{1}.intl_con,model_flux2{1}.biomass,'v:', 'Color',color_2m2, 'LineWidth',lineWidth); hold off
    xlim([min(model_flux2{1}.intl_con) max(model_flux2{1}.intl_con)])
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux'], 'FontSize',titleSize)
    legend({'2 Models, Model #1','2 Models, Model #2'}, ...
        'FontSize',legendSize, 'Location','SouthEast', 'EdgeColor','w');
    
    % Only Intracellular Sparsity Constraints With Growth and No Plateau
    model2_max = model_flux2{1}.intl_con(intlCon_idx2(find(abs(diff(model_flux2{1}.biomass(intlCon_idx2))) <= 1E-9,1))); if isempty(model2_max); model2_max = model_flux2{1}.intl_con(end); end
    model2_min = model_flux2{1}.intl_con(find(model_flux2{1}.biomass == 0, 1, 'last')); if isempty(model2_min); model2_min = model_flux2{1}.intl_con(1); end
    fig_handles(2) = figure(2);
    plot(model_flux2{1}.intl_con,model_flux2{1}.biomass,'^--', 'Color',color_2m1, 'LineWidth',lineWidth); hold on
    plot(model_flux2{1}.intl_con,model_flux2{1}.biomass,'v:', 'Color',color_2m2, 'LineWidth',lineWidth); hold off
    xlim([model2_min-1 model2_max+1])
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux'], 'FontSize',titleSize)
    legend({'2 Models, Model #1','2 Models, Model #2'}, ...
        'FontSize',legendSize, 'Location','SouthEast', 'EdgeColor','w');
    
elseif plotFlag == 1 % 1 Model
    % All Intracellular Sparsity Constraints
    fig_handles(1) = figure(1);
    plot(model_flux1{1}.intl_con,model_flux1{1}.biomass,'s-', 'Color',color_1m, 'LineWidth',lineWidth);
    xlim([min(model_flux1{1}.intl_con) max(model_flux1{1}.intl_con)])
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux'], 'FontSize',titleSize)
    legend({'1 Model'}, 'FontSize',legendSize, 'Location','SouthEast', 'EdgeColor','w');
    
    % Only Intracellular Sparsity Constraints With Growth and No Plateau
    model1_max = model_flux1{1}.intl_con(intlCon_idx1(find(abs(diff(model_flux1{1}.biomass(intlCon_idx1))) <= 1E-9,1))); if isempty(model1_max); model1_max = model_flux1{1}.intl_con(end); end
    model1_min = model_flux1{1}.intl_con(find(model_flux1{1}.biomass == 0, 1, 'last')); if isempty(model1_min); model1_min = model_flux1{1}.intl_con(1); end
    fig_handles(2) = figure(2);
    plot(model_flux1{1}.intl_con,model_flux1{1}.biomass,'s-', 'Color',color_1m, 'LineWidth',lineWidth);
    xlim([model1_min-1 model1_max+1])
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux'], 'FontSize',titleSize)
    legend({'1 Model'}, 'FontSize',legendSize, 'Location','SouthEast', 'EdgeColor','w');
end

%% Flux (Intracellular)
load('flux_color.mat')
tol = 1E-6;

if plotFlag == 3 % 1 Model & 2 Models
    % 1 Model
    data_flux_1m = model_flux1{1}.flux(model_flux1{1}.intl_idx,intlCon_idx1);
    data_int_1m = model_flux1{1}.int(model_flux1{1}.intl_idx,intlCon_idx1);
    % 2 Models, Model 1
    data_flux_2m1 = model_flux2{1}.flux(model_flux2{1}.intl_idx,intlCon_idx2);
    data_int_2m1 = model_flux2{1}.int(model_flux2{1}.intl_idx,intlCon_idx2);
    % 2 Models, Model 2
    data_flux_2m2 = model_flux2{2}.flux(model_flux2{2}.intl_idx,intlCon_idx2);
    data_int_2m2 = model_flux2{2}.int(model_flux2{2}.intl_idx,intlCon_idx2);
    
    [rxnsIdx,~] = biclusterFlux([data_int_1m, data_int_2m1, data_int_2m2],'single','jaccard');
    flux_lim = max(abs([data_flux_1m(:); data_flux_2m1(:); data_flux_2m2(:)]));
    
elseif plotFlag == 2 % 2 Models
    % 2 Models, Model 1
    data_flux_2m1 = model_flux2{1}.flux(model_flux2{1}.intl_idx,intlCon_idx2);
    data_int_2m1 = model_flux2{1}.int(model_flux2{1}.intl_idx,intlCon_idx2);
    data_int_2m1 = zeros(size(data_flux_2m1)); data_int_2m1(data_flux_2m1 ~= 0) = 1;
    % 2 Models, Model 2
    data_flux_2m2 = model_flux2{2}.flux(model_flux2{2}.intl_idx,intlCon_idx2);
    data_int_2m2 = model_flux2{2}.int(model_flux2{2}.intl_idx,intlCon_idx2);
    data_int_2m2 = zeros(size(data_flux_2m2)); data_int_2m2(data_flux_2m2 ~= 0) = 1;
    
    [rxnsIdx,~] = biclusterFlux([data_int_2m1, data_int_2m2],'single','jaccard');
    flux_lim = max(abs([data_flux_2m1(:); data_flux_2m2(:)]));
    
elseif plotFlag == 1 % 1 Model
    % 1 Model
    data_flux_1m = model_flux1{1}.flux(model_flux1{1}.intl_idx,intlCon_idx1);
    data_int_1m = model_flux1{1}.int(model_flux1{1}.intl_idx,intlCon_idx1);
    data_int_1m = zeros(size(data_flux_1m)); data_int_1m(abs(data_flux_1m) <= tol) = 1;

    [rxnsIdx,~] = biclusterFlux(data_int_1m,'single','jaccard');
    flux_lim = max(abs(data_flux_1m(:)));
end

if plotFlag == 3 || plotFlag == 1
    % 1 Model
    fig_handles(3) = figure(3);
    imagesc(data_flux_1m(rxnsIdx,:)); colormap(flux_color); caxis(flux_lim.*[-1 1]);
    cbar = colorbar; cbar.FontSize = 0.75*axesLabelSize;
    cbar.Label.String = 'Flux (mmol/gCDW/hr)'; cbar.Label.FontSize = xyLabelSize;
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45; ax.YTick = [];
    ax.XTick = 1:numel(intlCon_idx1); ax.XTickLabel = model_flux1{1}.intl_con(intlCon_idx1);
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Intracellular Reactions', 'FontSize',xyLabelSize)
    title([model_title ': 1 Model, Intracellular Flux'], 'FontSize',titleSize)
    
end

if plotFlag == 3 || plotFlag == 2
    % 2 Models, Model 1
    fig_handles(4) = figure(4);
    imagesc(data_flux_2m1(rxnsIdx,:)); colormap(flux_color); caxis(flux_lim.*[-1 1]);
    cbar = colorbar; cbar.FontSize = 0.75*axesLabelSize;
    cbar.Label.String = 'Flux (mmol/gCDW/hr)'; cbar.Label.FontSize = xyLabelSize;
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45; ax.YTick = [];
    ax.XTick = 1:numel(intlCon_idx2); ax.XTickLabel = model_flux2{1}.intl_con(intlCon_idx2);
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Intracellular Reactions', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Model 1 Intracellular Flux'], 'FontSize',titleSize)
    
    % 2 Models, Model 2
    fig_handles(5) = figure(5);
    imagesc(data_flux_2m2(rxnsIdx,:)); colormap(flux_color); caxis(flux_lim.*[-1 1]);
    cbar = colorbar; cbar.FontSize = 0.75*axesLabelSize;
    cbar.Label.String = 'Flux (mmol/gCDW/hr)'; cbar.Label.FontSize = xyLabelSize;
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45; ax.YTick = [];
    ax.XTick = 1:numel(intlCon_idx2); ax.XTickLabel = model_flux2{2}.intl_con(intlCon_idx2);
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Intracellular Reactions', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Model 2 Intracellular Flux'], 'FontSize',titleSize)
    
end

%% Binary Variables (Intracellular)

if plotFlag == 3 || plotFlag == 1
    % 1 Model
    fig_handles(6) = figure(6);
    imagesc(data_int_1m(rxnsIdx,:)); colormap(flip(gray(2))); caxis([0 1]);
    cbar = colorbar; cbar.Label.String = 'Binary Variable'; cbar.Label.FontSize = xyLabelSize;
    cbar.Ticks = [0.25 0.75]; cbar.TickLabels = [0,1]; cbar.FontSize = 0.75*axesLabelSize;
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45; ax.YTick = [];
    ax.XTick = 1:numel(intlCon_idx1); ax.XTickLabel = model_flux1{1}.intl_con(intlCon_idx1);
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Intracellular Reactions', 'FontSize',xyLabelSize)
    title([model_title ': 1 Model, Intracellular Binary Variables'], 'FontSize',titleSize)
end

if plotFlag == 3 || plotFlag == 2
    % 2 Models, Model 1
    fig_handles(7) = figure(7);
    imagesc(data_int_2m1(rxnsIdx,:)); colormap(flip(gray(2))); caxis([0 1]);
    cbar = colorbar; cbar.Label.String = 'Binary Variable'; cbar.Label.FontSize = xyLabelSize;
    cbar.Ticks = [0.25 0.75]; cbar.TickLabels = [0,1]; cbar.FontSize = 0.75*axesLabelSize;
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45; ax.YTick = [];
    ax.XTick = 1:numel(intlCon_idx2); ax.XTickLabel = model_flux2{1}.intl_con(intlCon_idx2);
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Intracellular Reactions', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Model 1 Intracellular Binary Variables'], 'FontSize',titleSize)
    
    % 2 Models, Model 2
    fig_handles(8) = figure(8);
    imagesc(data_int_2m2(rxnsIdx,:)); colormap(flip(gray(2))); caxis([0 1]);
    cbar = colorbar; cbar.Label.String = 'Binary Variable'; cbar.Label.FontSize = xyLabelSize;
    cbar.Ticks = [0.25 0.75]; cbar.TickLabels = [0,1]; cbar.FontSize = 0.75*axesLabelSize;
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45; ax.YTick = [];
    ax.XTick = 1:numel(intlCon_idx2); ax.XTickLabel = model_flux2{2}.intl_con(intlCon_idx2);
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Intracellular Reactions', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Model 2 Intracellular Binary Variables'], 'FontSize',titleSize)
end

%% Dendrograms
col_thresh = 0.33;

if plotFlag == 3 || plotFlag == 1
    % 1 Model
    Z = linkage(single(data_int_1m)','single','jaccard');
    fig_handles(9) = figure(9);
    dendrogram(Z,0, 'labels',int2str(model_flux1{1}.intl_con(intlCon_idx1)'), 'Reorder',1:numel(intlCon_idx1), 'ColorThreshold',col_thresh*max(Z(:,3)), 'CheckCrossing',false); hold on
    plot(0.5+[0 numel(intlCon_idx1)],col_thresh*max(Z(:,3)).*ones(2,1),'k--', 'LineWidth',0.5*lineWidth); hold off
    xlim(0.5+[0 numel(intlCon_idx1)]); ylim([0 1.1*max(Z(:,3))]); box on;
    ax = gca; ax.FontSize = 0.75*axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Jaccard Distance', 'FontSize',xyLabelSize)
    title([model_title ': 1 Model'], 'FontSize',titleSize)
end

if plotFlag == 3 || plotFlag == 2
    % 2 Models, Model 1
    Z = linkage(single(data_int_2m1)','single','jaccard');
    fig_handles(10) = figure(10);
    dendrogram(Z,0, 'Labels',int2str(model_flux2{1}.intl_con(intlCon_idx2)'), 'Reorder',1:numel(intlCon_idx2), 'ColorThreshold',col_thresh*max(Z(:,3)), 'CheckCrossing',false); hold on
    plot(0.5+[0 numel(intlCon_idx2)],col_thresh*max(Z(:,3)).*ones(2,1),'k--', 'LineWidth',0.5*lineWidth); hold off
    xlim(0.5+[0 numel(intlCon_idx2)]); ylim([0 1.1*max(Z(:,3))]); box on;
    ax = gca; ax.FontSize = 0.75*axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Jaccard Distance', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Model 1'], 'FontSize',titleSize)
    
    % 2 Models, Model 2
    Z = linkage(single(data_int_2m2)','single','jaccard');
    fig_handles(11) = figure(11);
    dendrogram(Z,0, 'Labels',int2str(model_flux2{2}.intl_con(intlCon_idx2)'), 'Reorder',1:numel(intlCon_idx2), 'ColorThreshold',col_thresh*max(Z(:,3)), 'CheckCrossing',false); hold on
    plot(0.5+[0 numel(intlCon_idx2)],col_thresh*max(Z(:,3)).*ones(2,1),'k--', 'LineWidth',0.5*lineWidth); hold off
    xlim(0.5+[0 numel(intlCon_idx2)]); ylim([0 1.1*max(Z(:,3))]); box on;
    ax = gca; ax.FontSize = 0.75*axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Jaccard Distance', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Model 2'], 'FontSize',titleSize)
end

faceAlpha = 0.1;
if plotFlag == 3
    % 1 Model & 2 Models
    data_int = [data_int_1m, data_int_2m1, data_int_2m2];
    labels = [arrayfun(@(x) [int2str(model_flux1{1}.intl_con(intlCon_idx1(x))) '_1m'],1:numel(intlCon_idx1), 'Uni',false),...
        arrayfun(@(x) [int2str(model_flux2{1}.intl_con(intlCon_idx2(x))) '_2m1'],1:numel(intlCon_idx2), 'Uni',false), ...
        arrayfun(@(x) [int2str(model_flux2{2}.intl_con(intlCon_idx2(x))) '_2m2'],1:numel(intlCon_idx2), 'Uni',false)];
    Z = linkage(single(data_int)','single','jaccard');
    fig_handles(12) = figure(12);
    temp = Z(:,3); temp(temp == max(temp)) = [];
    dendrogram(Z,0, 'Labels',labels, 'Reorder',1:numel(labels), 'ColorThreshold',col_thresh*max(temp), 'CheckCrossing',false); hold on
    plot(0.5+[0 numel(intlCon_idx1)+numel(intlCon_idx2)*2],col_thresh*max(temp).*ones(2,1),'k--', 'LineWidth',0.5*lineWidth);
    area(0.5+[0, numel(intlCon_idx1)],1.1*max(temp).*[1, 1], 'FaceColor',color_1m, 'EdgeAlpha',0, 'FaceAlpha',faceAlpha);
    area(0.5+[numel(intlCon_idx1), numel(intlCon_idx1)+numel(intlCon_idx2)],1.1*max(temp).*[1, 1], 'FaceColor',color_2m1, 'EdgeAlpha',0, 'FaceAlpha',faceAlpha);
    area(0.5+[numel(intlCon_idx1)+numel(intlCon_idx2), numel(intlCon_idx1)+numel(intlCon_idx2)*2],1.1*max(temp).*[1, 1], 'FaceColor',color_2m2, 'EdgeAlpha',0, 'FaceAlpha',faceAlpha); hold off
    xlim(0.5+[0 numel(intlCon_idx1)+numel(intlCon_idx2)*2]); ylim([0 1.1*max(temp)]); box on;
    ax = gca; ax.FontSize = 0.75*axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Jaccard Distance', 'FontSize',xyLabelSize)
    title([model_title ': 1 Model & 2 Models'], 'FontSize',titleSize)
    
elseif plotFlag == 2
    % 2 Models
    data_int = [data_int_2m1, data_int_2m2];
    labels = [arrayfun(@(x) [int2str(model_flux2{1}.intl_con(intlCon_idx2(x))) '_2m1'],1:numel(intlCon_idx2), 'Uni',false), ...
        arrayfun(@(x) [int2str(model_flux2{2}.intl_con(intlCon_idx2(x))) '_2m2'],1:numel(intlCon_idx2), 'Uni',false)];
    Z = linkage(single(data_int)','single','jaccard');
    fig_handles(12) = figure(12);
    temp = Z(:,3); temp(temp == max(temp)) = [];
    dendrogram(Z,0, 'Labels',labels, 'Reorder',1:numel(labels), 'ColorThreshold',col_thresh*max(temp), 'CheckCrossing',false); hold on
    plot(0.5+[0 numel(intlCon_idx2)*2],col_thresh*max(temp).*ones(2,1),'k--', 'LineWidth',0.5*lineWidth);
    area(0.5+[0, numel(intlCon_idx2)],1.1*max(temp).*[1, 1], 'FaceColor',color_2m1, 'EdgeAlpha',0, 'FaceAlpha',faceAlpha);
    area(0.5+[numel(intlCon_idx2), numel(intlCon_idx2)*2],1.1*max(temp).*[1, 1], 'FaceColor',color_2m2, 'EdgeAlpha',0, 'FaceAlpha',faceAlpha); hold off
    xlim(0.5+[0 numel(intlCon_idx2)*2]); ylim([0 1.1*max(temp)]); box on;
    ax = gca; ax.FontSize = 0.75*axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Jaccard Distance', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Clustered by Intracellular Flux'], 'FontSize',titleSize)
end

%% T-SNE

if plotFlag == 3
    data_flux = [data_flux_1m, data_flux_2m1, data_flux_2m2];
    idx_1m = 1:numel(intlCon_idx1);
    idx_2m1 = numel(intlCon_idx1)+(1:numel(intlCon_idx2));
    idx_2m2 = (numel(intlCon_idx1)+numel(intlCon_idx2))+(1:numel(intlCon_idx2));
    
    options = statset('MaxIter',1E4);
    tsne_loss = 1; iter = 1;
    while tsne_loss > 0.25 && iter < 100
        [tsne_flux,tsne_loss] = tsne(data_flux', 'Algorithm','exact', 'Options',options, 'Distance','spearman', 'Perplexity',min([numel(intlCon_idx1), numel(intlCon_idx2)]));
        iter = iter + 1;
    end
    intlCon = union(model_flux1{1}.intl_con(intlCon_idx1),model_flux2{1}.intl_con(intlCon_idx2));
    
    tsne_color = flip(gray(numel(intlCon)-1));
    fig_handles(13) = figure(13);
    h(1) = scatter(tsne_flux(idx_1m(end),1),tsne_flux(idx_1m(end),2),markerSize,color_1m,'s','filled', 'LineWidth',lineWidth, 'MarkerEdgeColor','k'); hold on
    h(2) = scatter(tsne_flux(idx_2m1(end),1),tsne_flux(idx_2m1(end),2),0.75*markerSize,color_2m1,'^','filled', 'LineWidth',lineWidth, 'MarkerEdgeColor','k');
    h(3) = scatter(tsne_flux(idx_2m2(end),1),tsne_flux(idx_2m2(end),2),0.75*markerSize,color_2m2,'v','filled', 'LineWidth',lineWidth, 'MarkerEdgeColor','k');
    for ii = 1:numel(intlCon)-1
        [~,~,idx] = intersect(intlCon(ii),model_flux1{1}.intl_con(intlCon_idx1));
        if ~isempty(idx)
            scatter(tsne_flux(idx_1m(idx),1),tsne_flux(idx_1m(idx),2),markerSize,tsne_color(ii,:),'s','filled', 'LineWidth',0.5*lineWidth, 'MarkerEdgeColor','k');
        end
        [~,~,idx] = intersect(intlCon(ii),model_flux2{1}.intl_con(intlCon_idx2));
        if ~isempty(idx)
            scatter(tsne_flux(idx_2m1(idx),1),tsne_flux(idx_2m1(idx),2),0.75*markerSize,tsne_color(ii,:),'^','filled', 'LineWidth',0.5*lineWidth, 'MarkerEdgeColor','k');
            scatter(tsne_flux(idx_2m2(idx),1),tsne_flux(idx_2m2(idx),2),0.75*markerSize,tsne_color(ii,:),'v','filled', 'LineWidth',0.5*lineWidth, 'MarkerEdgeColor','k');
        end
    end
    hold off; axis tight; ax = gca; ax.Visible = 'off';
    colormap(tsne_color); caxis([1, numel(intlCon)]);
    cbar = colorbar; cbar.FontSize = 0.75*axesLabelSize; cbar.TickDirection = 'out';
    cbar.Ticks = 0.5+(1:numel(intlCon)-1); cbar.TickLabels = intlCon(1:end-1);
    cbar.Label.String = 'Intracellular Sparsity Constraint'; cbar.Label.FontSize = xyLabelSize;
    legend(h,{'1 Model','2 Models, Model #1','2 Models, Model #2'}, 'FontSize',legendSize, 'Location','Best', 'EdgeColor','w');
    title([model_title ': 1 Model & 2 Models, Intracellular Flux t-SNE'], 'FontSize',titleSize)
    set(findall(ax,'type','text'), 'visible','on')
elseif plotFlag == 2
    data_flux = [data_flux_2m1, data_flux_2m2];
    idx_2m1 = 1:numel(intlCon_idx2);
    idx_2m2 = numel(intlCon_idx2)+(1:numel(intlCon_idx2));
    
    options = statset('MaxIter',1E4);
    tsne_loss = 1; iter = 1;
    while tsne_loss > 0.25 && iter < 100
        [tsne_flux,tsne_loss] = tsne(data_flux', 'Algorithm','exact', 'Options',options, 'Distance','spearman', 'Perplexity',numel(intlCon_idx2));
        iter = iter + 1;
    end
    intlCon = model_flux2{1}.intl_con(intlCon_idx2);
    
    tsne_color = flip(gray(numel(intlCon)-1));
    fig_handles(13) = figure(13);
    h(1) = scatter(tsne_flux(idx_2m1(end),1),tsne_flux(idx_2m1(end),2),0.75*markerSize,color_2m1,'^','filled', 'LineWidth',lineWidth, 'MarkerEdgeColor','k'); hold on
    h(2) = scatter(tsne_flux(idx_2m2(end),1),tsne_flux(idx_2m2(end),2),0.75*markerSize,color_2m2,'v','filled', 'LineWidth',lineWidth, 'MarkerEdgeColor','k');
    for ii = 1:numel(intlCon)-1
        [~,~,idx] = intersect(intlCon(ii),model_flux2{1}.intl_con(intlCon_idx2));
        if ~isempty(idx)
            scatter(tsne_flux(idx_2m1(idx),1),tsne_flux(idx_2m1(idx),2),0.75*markerSize,tsne_color(ii,:),'^','filled', 'LineWidth',0.5*lineWidth, 'MarkerEdgeColor','k');
            scatter(tsne_flux(idx_2m2(idx),1),tsne_flux(idx_2m2(idx),2),0.75*markerSize,tsne_color(ii,:),'v','filled', 'LineWidth',0.5*lineWidth, 'MarkerEdgeColor','k');
        end
    end
    hold off; axis tight; ax = gca; ax.Visible = 'off';
    colormap(tsne_color); caxis([1, numel(intlCon)]);
    cbar = colorbar; cbar.FontSize = 0.75*axesLabelSize; cbar.TickDirection = 'out';
    cbar.Ticks = 0.5+(1:numel(intlCon)-1); cbar.TickLabels = intlCon(1:end-1);
    cbar.Label.String = 'Intracellular Sparsity Constraint'; cbar.Label.FontSize = xyLabelSize;
    legend(h,{'2 Models, Model #1','2 Models, Model #2'}, 'FontSize',legendSize, 'Location','Best', 'EdgeColor','w');
    title([model_title ': 2 Models, Intracellular Flux t-SNE'], 'FontSize',titleSize)
    set(findall(ax,'type','text'), 'visible','on')
end

%% Flux Correlations (Intracellular)

if plotFlag == 3 || plotFlag == 2
    R = corr(data_flux_2m1',data_flux_2m2', 'Type','Spearman'); R(isnan(R)) = 0;
    [P1,P2] = biclusterFlux(R); R = R(P1,P2);
    rxns_2m1 = model_flux2{1}.rxns(model_flux2{1}.intl_idx(P1));
    idx_2m1 = find(ismember(R,zeros(1,size(R,1)),'rows'));
    R(idx_2m1,:) = []; rxns_2m1(idx_2m1) = [];
    rxns_2m2 = model_flux2{2}.rxns(model_flux2{2}.intl_idx(P2));
    idx_2m2 = find(ismember(R',zeros(1,size(R',2)),'rows'));
    R(:,idx_2m2) = []; rxns_2m2(idx_2m2) = [];
    
    fig_handles(14) = figure(14);
    imagesc(R'); colormap(flux_color); caxis([-1 1]); axis square; grid off
    cbar = colorbar; cbar.FontSize = 0.75*axesLabelSize;
    cbar.Label.String = 'Spearman Correlation'; cbar.Label.FontSize = xyLabelSize;
    ax = gca; ax.FontSize = 0.25*axesLabelSize; ax.XTickLabelRotation = 90;
    ax.XTick = 1:numel(rxns_2m1); ax.XTickLabel = rxns_2m1;
    ax.YTick = 1:numel(rxns_2m2); ax.YTickLabel = rxns_2m2;
    ax.TickLength = tickLength.*[0.5 1];
    xlabel('Model 1 Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Model 2 Intracellular Reactions', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Intracellular Flux'], 'FontSize',titleSize)

end

%% Flux (Exchange)

if plotFlag == 3
    % 1 Model
    data_flux_1m = model_flux1{1}.flux(model_flux1{1}.exch_idx,intlCon_idx1);
    % 2 Models, Model 1
    data_flux_2m1 = model_flux2{1}.flux(model_flux2{1}.exch_idx,intlCon_idx2);
    % 2 Models, Model 2
    data_flux_2m2 = model_flux2{2}.flux(model_flux2{2}.exch_idx,intlCon_idx2);
    
    [rxnsIdx,~] = biclusterFlux([data_flux_1m, data_flux_2m1, data_flux_2m2]);
    rm_idx = find(ismember([data_flux_1m(rxnsIdx,:), data_flux_2m1(rxnsIdx,:), data_flux_2m2(rxnsIdx,:)],...
        zeros(1,numel([intlCon_idx1, intlCon_idx2, intlCon_idx2])),'rows'));
    rxnsIdx(rm_idx) = []; clear rm_idx
    rxns = strrep(model_flux1{1}.rxnNames(model_flux1{1}.exch_idx(rxnsIdx)),' exchange','');
    
    flux_lim = max(abs([data_flux_1m(:); data_flux_2m1(:); data_flux_2m2(:)]));
    
elseif plotFlag == 2
    % 2 Models, Model 1
    data_flux_2m1 = model_flux2{1}.flux(model_flux2{1}.exch_idx,intlCon_idx2);
    % 2 Models, Model 2
    data_flux_2m2 = model_flux2{2}.flux(model_flux2{2}.exch_idx,intlCon_idx2);
    
    [rxnsIdx,~] = biclusterFlux([data_flux_2m1, data_flux_2m2]);
    rm_idx = find(ismember([data_flux_2m1(rxnsIdx,:), data_flux_2m2(rxnsIdx,:)],...
        zeros(1,numel([intlCon_idx2, intlCon_idx2])),'rows'));
    rxnsIdx(rm_idx) = []; clear rm_idx
    rxns = strrep(model_flux2{1}.rxnNames(model_flux2{1}.exch_idx(rxnsIdx)),' exchange','');
    
    flux_lim = max(abs([data_flux_2m1(:); data_flux_2m2(:)]));
    
elseif plotFlag == 1
    % 1 Model
    data_flux_1m = model_flux1{1}.flux(model_flux1{1}.exch_idx,intlCon_idx1);
    
    [rxnsIdx,~] = biclusterFlux(data_flux_1m);
    rm_idx = find(ismember(data_flux_1m(rxnsIdx,:),...
        zeros(1,numel(intlCon_idx1)),'rows'));
    rxnsIdx(rm_idx) = []; clear rm_idx
    rxns = strrep(model_flux1{1}.rxnNames(model_flux1{1}.exch_idx(rxnsIdx)),' exchange','');
    
    flux_lim = max(abs(data_flux_1m(:)));
end

if plotFlag == 3 || plotFlag == 1
    % 1 Model
    fig_handles(15) = figure(15);
    imagesc(data_flux_1m(rxnsIdx,:)); colormap(flux_color); caxis(flux_lim.*[-1 1]);
    cbar = colorbar; cbar.FontSize = 0.75*axesLabelSize;
    cbar.Label.String = 'Flux (mmol/gCDW/hr)'; cbar.Label.FontSize = xyLabelSize;
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    ax.XTick = 1:numel(intlCon_idx1); ax.XTickLabel = model_flux1{1}.intl_con(intlCon_idx1);
    ax.YTick = 1:numel(rxns); ax.YTickLabel = rxns;
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Exchange Reactions', 'FontSize',xyLabelSize)
    title([model_title ': 1 Model, Exchange Flux'], 'FontSize',titleSize)
    
end

if plotFlag == 3 || plotFlag == 2
    % 2 Models, Model 1
    fig_handles(16) = figure(16);
    imagesc(data_flux_2m1(rxnsIdx,:)); colormap(flux_color); caxis(flux_lim.*[-1 1]);
    cbar = colorbar; cbar.FontSize = 0.75*axesLabelSize;
    cbar.Label.String = 'Flux (mmol/gCDW/hr)'; cbar.Label.FontSize = xyLabelSize;
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    ax.XTick = 1:numel(intlCon_idx2); ax.XTickLabel = model_flux2{1}.intl_con(intlCon_idx2);
    ax.YTick = 1:numel(rxns); ax.YTickLabel = rxns;
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Exchange Reactions', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Model 1 Exchange Flux'], 'FontSize',titleSize)
    
    % 2 Models, Model 2
    fig_handles(17) = figure(17);
    imagesc(data_flux_2m2(rxnsIdx,:)); colormap(flux_color); caxis(flux_lim.*[-1 1]);
    cbar = colorbar; cbar.FontSize = 0.75*axesLabelSize;
    cbar.Label.String = 'Flux (mmol/gCDW/hr)'; cbar.Label.FontSize = xyLabelSize;
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    ax.XTick = 1:numel(intlCon_idx2); ax.XTickLabel = model_flux2{2}.intl_con(intlCon_idx2);
    ax.YTick = 1:numel(rxns); ax.YTickLabel = rxns;
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Exchange Reactions', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Model 2 Exchange Flux'], 'FontSize',titleSize)

end

%% Flux Correlations (Exchange)

if plotFlag == 3 || plotFlag == 2
    R = corr(data_flux_2m1',data_flux_2m2', 'Type','Spearman'); R(isnan(R)) = 0;
    [P1,P2] = biclusterFlux(R); R = R(P1,P2);
    rxns_2m1 = strrep(model_flux2{1}.rxnNames(model_flux2{1}.exch_idx(P1)),' exchange','');
    idx_2m1 = find(ismember(R,zeros(1,size(R,1)),'rows'));
    R(idx_2m1,:) = []; rxns_2m1(idx_2m1) = [];
    rxns_2m2 = strrep(model_flux2{2}.rxnNames(model_flux2{2}.exch_idx(P2)),' exchange','');
    idx_2m2 = find(ismember(R',zeros(1,size(R',2)),'rows'));
    R(:,idx_2m2) = []; rxns_2m2(idx_2m2) = [];
    
    fig_handles(18) = figure(18);
    imagesc(R'); colormap(flux_color); caxis([-1 1]); axis square; grid off
    cbar = colorbar; cbar.FontSize = 0.75*axesLabelSize;
    cbar.Label.String = 'Spearman Correlation'; cbar.Label.FontSize = xyLabelSize;
    ax = gca; ax.FontSize = 0.75*axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.XTickLabelRotation = 45; ax.TickLabelInterpreter = 'none';
    ax.XTick = 1:numel(rxns_2m1); ax.XTickLabel = rxns_2m1;
    ax.YTick = 1:numel(rxns_2m2); ax.YTickLabel = rxns_2m2;
    xlabel('Model 1 Exchange Reactions', 'FontSize',xyLabelSize)
    ylabel('Model 2 Exchange Reactions', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Exchange Flux'], 'FontSize',titleSize)
end

%% Extracellular Metabolites

if plotFlag == 3 || plotFlag == 2
    load('classifyExchangeFlux_cmap.mat')
    
    % Metabolites that are Taken Up, Produced, or Exchanged
    [exchMets_idx,~] = find(model_flux2{1}.exchMets.fluxType(:,intlCon_idx2) ~= 0); exchMets_idx = unique(exchMets_idx);
    exchMets = erase(model_flux2{1}.rxnNames(model_flux2{1}.exch_idx(exchMets_idx)),' exchange');
    fig = figure(19); ax = gca;
    [fig,ax,cbar] = plotClassifyExchangeFlux(model_flux2{1}.exchMets.fluxType(exchMets_idx,intlCon_idx2),fig,ax,exchMets);
    ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    ax.XTick = 1:numel(intlCon_idx2); ax.XTickLabel = model_flux2{1}.intl_con(intlCon_idx2);
    cbar.FontSize = 0.75*axesLabelSize; cbar.Label.FontSize = xyLabelSize;
    xlabel(ax, 'Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel(ax, 'Extracellular Metabolites', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Extracellular Metabolites with Non-Zero Flux'], 'FontSize',titleSize)
    box on; grid on
    fig_handles(19) = fig;
    
    % Metabolites that are Exchanged
    [exchMets_idx,~] = find(model_flux2{1}.exchMets.model1_to_model2 + model_flux2{1}.exchMets.model2_to_model1); exchMets_idx = unique(exchMets_idx);
    exchMets = erase(model_flux2{1}.rxnNames(model_flux2{1}.exch_idx(exchMets_idx)),' exchange');
    fig = figure(20); ax = gca;
    [fig,ax,cbar] = plotClassifyExchangeFlux(model_flux2{1}.exchMets.fluxType(exchMets_idx,intlCon_idx2),fig,ax,exchMets);
    ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    ax.XTick = 1:numel(intlCon_idx2); ax.XTickLabel = model_flux2{1}.intl_con(intlCon_idx2);
    cbar.FontSize = 0.75*axesLabelSize; cbar.Label.FontSize = xyLabelSize;
    xlabel(ax, 'Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel(ax, 'Extracellular Metabolites', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Extracellular Metabolites that are Exchanged'], 'FontSize',titleSize)
    box on; grid on
    fig_handles(20) = fig;
    
    % Number of Times Metabolites are Exchanged
    fig_handles(21) = figure(21);
    [~,sorted_idx] = sort(model_flux2{1}.exchMets.model1_to_model2(exchMets_idx,:)+model_flux2{1}.exchMets.model2_to_model1(exchMets_idx,:));
    h = barh(100.*[sum(model_flux2{1}.exchMets.model1_to_model2(exchMets_idx,:),2), ...
        sum(model_flux2{1}.exchMets.model2_to_model1(exchMets_idx,:),2)]./numel(intlCon_idx2),'stacked');
    h(1).FaceColor = fluxType_cmap(8,:); h(2).FaceColor = fluxType_cmap(9,:);
    xlim([0 100]); ylim([0 numel(exchMets)+1])
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none';
    ax.YTick = 1:numel(exchMets); ax.YTickLabels = exchMets(sorted_idx);
    legend({'2\rightarrow1','1\rightarrow2'}, 'FontSize',legendSize, 'Location','SouthEast', 'EdgeColor','w')
    xlabel('Percentage of Times Exchanged', 'FontSize',xyLabelSize)
    ylabel('Extracellular Metabolite', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Extracellular Metabolites that are Exchanged'], 'FontSize',titleSize)
    
    % Number of Metabolites Exchanged
    fig_handles(22) = figure(22);
    num_12 = sum(model_flux2{1}.exchMets.model1_to_model2)';
    num_21 = sum(model_flux2{1}.exchMets.model2_to_model1)';
    h = bar([num_21(intlCon_idx2), num_12(intlCon_idx2)],'stacked');
    h(1).FaceColor = fluxType_cmap(8,:); h(2).FaceColor = fluxType_cmap(9,:);
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    ax.XTick = 1:numel(intlCon_idx2); ax.XTickLabel = model_flux2{1}.intl_con(intlCon_idx2);
    ytick_idx = (ax.YTick == round(ax.YTick)); ax.YTick = ax.YTick(ytick_idx);
    xlim([0 numel(intlCon_idx2)+1]); ylim([0 max(num_21(intlCon_idx2) + num_12(intlCon_idx2))+0.5])
    legend({'2\rightarrow1','1\rightarrow2'}, 'FontSize',legendSize, 'Location','NorthEast', 'EdgeColor','w')
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Number of Metabolites Exchanged', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Extracellular Metabolites that are Exchanged'], 'FontSize',titleSize)
end

%% FBA Solutions

if plotFlag == 3
    diff_1m = model_flux1{1}.warningFlag.sol - model_flux1{1}.biomass;
    diff_2m1 = model_flux2{1}.warningFlag.sol - model_flux2{1}.biomass;
    diff_2m2 = model_flux2{2}.warningFlag.sol - model_flux2{2}.biomass;
    
    max_diff = max([max(abs(diff_1m(:))), max(abs(diff_2m1(:))), max(abs(diff_2m2(:)))]);
    
    % All Intracellular Sparsity Constraints
    fig_handles(23) = figure(23);
    plot(model_flux1{1}.intl_con,diff_1m,'s-', 'Color',color_1m, 'LineWidth',lineWidth); hold on
    plot(model_flux2{1}.intl_con,diff_2m1,'^--', 'Color',color_2m1, 'LineWidth',lineWidth);
    plot(model_flux2{1}.intl_con,diff_2m2,'v:', 'Color',color_2m2, 'LineWidth',lineWidth); hold off
    if sum(model_flux1{1}.warningFlag.FBA_sol) % if model cannot grow, mark with a black x
        hold on; plot(model_flux1{1}.intl_con(find(model_flux1{1}.warningFlag.FBA_sol)),diff_1m(find(model_flux1{1}.warningFlag.FBA_sol)),'kx', 'LineWidth',3); hold off
    end
    if sum(model_flux2{1}.warningFlag.FBA_sol) % if model cannot grow, mark with a black x
        hold on; plot(model_flux2{1}.intl_con(find(model_flux2{1}.warningFlag.FBA_sol)),diff_2m1(find(model_flux2{1}.warningFlag.FBA_sol)),'kx', 'LineWidth',3); hold off
    end
    if sum(model_flux2{2}.warningFlag.FBA_sol) % if model cannot grow, mark with a black x
        hold on; plot(model_flux2{2}.intl_con(find(model_flux2{2}.warningFlag.FBA_sol)),diff_2m2(find(model_flux2{2}.warningFlag.FBA_sol)),'kx', 'LineWidth',3); hold off
    end
    xlim([min([model_flux1{1}.intl_con, model_flux2{1}.intl_con]) max([model_flux1{1}.intl_con, model_flux2{1}.intl_con])])
    ylim(max_diff.*[-1 1])
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Difference in Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux (FBA - Algorithm)'], 'FontSize',titleSize)
    legend({'1 Model','2 Models, Model #1','2 Models, Model #2'}, 'FontSize',legendSize, 'Location','SouthEast', 'EdgeColor','w');
    
    % Only Intracellular Sparsity Constraints With Growth and No Plateau
    fig_handles(24) = figure(24);
    plot(model_flux1{1}.intl_con,diff_1m,'s-', 'Color',color_1m, 'LineWidth',lineWidth); hold on
    plot(model_flux2{1}.intl_con,diff_2m1,'^--', 'Color',color_2m1, 'LineWidth',lineWidth);
    plot(model_flux2{1}.intl_con,diff_2m2,'v:', 'Color',color_2m2, 'LineWidth',lineWidth); hold off
    if sum(model_flux1{1}.warningFlag.FBA_sol) % if model cannot grow, mark with a black x
        hold on; plot(model_flux1{1}.intl_con(find(model_flux1{1}.warningFlag.FBA_sol)),diff_1m(find(model_flux1{1}.warningFlag.FBA_sol)),'kx', 'LineWidth',3); hold off
    end
    if sum(model_flux2{1}.warningFlag.FBA_sol) % if model cannot grow, mark with a black x
        hold on; plot(model_flux2{1}.intl_con(find(model_flux2{1}.warningFlag.FBA_sol)),diff_2m1(find(model_flux2{1}.warningFlag.FBA_sol)),'kx', 'LineWidth',3); hold off
    end
    if sum(model_flux2{2}.warningFlag.FBA_sol) % if model cannot grow, mark with a black x
        hold on; plot(model_flux2{2}.intl_con(find(model_flux2{2}.warningFlag.FBA_sol)),diff_2m2(find(model_flux2{2}.warningFlag.FBA_sol)),'kx', 'LineWidth',3); hold off
    end
    xlim([nanmin(model1_min,model2_min)-1 nanmax(model1_max,model2_max)+1])
    ylim(max_diff.*[-1 1])
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Difference in Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux (FBA - Algorithm)'], 'FontSize',titleSize)
    legend({'1 Model','2 Models, Model #1','2 Models, Model #2'}, 'FontSize',legendSize, 'Location','SouthEast', 'EdgeColor','w');
    
elseif plotFlag == 2
    diff_2m1 = model_flux2{1}.warningFlag.sol - model_flux2{1}.biomass;
    diff_2m2 = model_flux2{2}.warningFlag.sol - model_flux2{2}.biomass;
    
    max_diff = max([max(abs(diff_2m1(:))), max(abs(diff_2m2(:)))]);
    
    % All Intracellular Sparsity Constraints
    fig_handles(23) = figure(23);
    plot(model_flux2{1}.intl_con,diff_2m1,'^--', 'Color',color_2m1, 'LineWidth',lineWidth); hold on
    plot(model_flux2{1}.intl_con,diff_2m2,'v:', 'Color',color_2m2, 'LineWidth',lineWidth); hold off
    if sum(model_flux2{1}.warningFlag.FBA_sol) % if model cannot grow, mark with a black x
        hold on; plot(model_flux2{1}.intl_con(find(model_flux2{1}.warningFlag.FBA_sol)),diff_2m1(find(model_flux2{1}.warningFlag.FBA_sol)),'kx', 'LineWidth',3); hold off
    end
    if sum(model_flux2{2}.warningFlag.FBA_sol) % if model cannot grow, mark with a black x
        hold on; plot(model_flux2{2}.intl_con(find(model_flux2{2}.warningFlag.FBA_sol)),diff_2m2(find(model_flux2{2}.warningFlag.FBA_sol)),'kx', 'LineWidth',3); hold off
    end
    xlim([min(model_flux2{1}.intl_con) max(model_flux2{1}.intl_con)])
    ylim(max_diff.*[-1 1])
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Difference in Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux (FBA - Algorithm)'], 'FontSize',titleSize)
    legend({'2 Models, Model #1','2 Models, Model #2'}, 'FontSize',legendSize, 'Location','SouthEast', 'EdgeColor','w');
    
    % Only Intracellular Sparsity Constraints With Growth and No Plateau
    fig_handles(24) = figure(24);
    plot(model_flux2{1}.intl_con,diff_2m1,'^--', 'Color',color_2m1, 'LineWidth',lineWidth); hold on
    plot(model_flux2{1}.intl_con,diff_2m2,'v:', 'Color',color_2m2, 'LineWidth',lineWidth); hold off
    if sum(model_flux2{1}.warningFlag.FBA_sol) % if model cannot grow, mark with a black x
        hold on; plot(model_flux2{1}.intl_con(find(model_flux2{1}.warningFlag.FBA_sol)),diff_2m1(find(model_flux2{1}.warningFlag.FBA_sol)),'kx', 'LineWidth',3); hold off
    end
    if sum(model_flux2{2}.warningFlag.FBA_sol) % if model cannot grow, mark with a black x
        hold on; plot(model_flux2{2}.intl_con(find(model_flux2{2}.warningFlag.FBA_sol)),diff_2m2(find(model_flux2{2}.warningFlag.FBA_sol)),'kx', 'LineWidth',3); hold off
    end
    xlim([nanmin(model2_min)-1 nanmax(model2_max)+1])
    ylim(max_diff.*[-1 1])
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Difference in Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux (FBA - Algorithm)'], 'FontSize',titleSize)
    legend({'2 Models, Model #1','2 Models, Model #2'}, 'FontSize',legendSize, 'Location','SouthEast', 'EdgeColor','w');
    
elseif plotFlag == 1
    diff_1m = model_flux1{1}.warningFlag.sol - model_flux1{1}.biomass;
    max_diff = max(abs(diff_1m(:)));
    
    % All Intracellular Sparsity Constraints
    fig_handles(23) = figure(23);
    plot(model_flux1{1}.intl_con,diff_1m,'s-', 'Color',color_1m, 'LineWidth',lineWidth);
    if sum(model_flux1{1}.warningFlag.FBA_sol) % if model cannot grow, mark with a black x
        hold on; plot(model_flux1{1}.intl_con(find(model_flux1{1}.warningFlag.FBA_sol)),diff_1m(find(model_flux1{1}.warningFlag.FBA_sol)),'kx', 'LineWidth',3); hold off
    end
    xlim([min(model_flux1{1}.intl_con) max(model_flux1{1}.intl_con)])
    ylim(max_diff.*[-1 1])
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Difference in Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux (FBA - Algorithm)'], 'FontSize',titleSize)
    legend({'1 Model','2 Models, Model #1','2 Models, Model #2'}, 'FontSize',legendSize, 'Location','SouthEast', 'EdgeColor','w');
    
    % Only Intracellular Sparsity Constraints With Growth and No Plateau
    fig_handles(24) = figure(24);
    plot(model_flux1{1}.intl_con,diff_1m,'s-', 'Color',color_1m, 'LineWidth',lineWidth);
    if sum(model_flux1{1}.warningFlag.FBA_sol) % if model cannot grow, mark with a black x
        hold on; plot(model_flux1{1}.intl_con(find(model_flux1{1}.warningFlag.FBA_sol)),diff_1m(find(model_flux1{1}.warningFlag.FBA_sol)),'kx', 'LineWidth',3); hold off
    end
    xlim([nanmin(model1_min)-1 nanmax(model1_max)+1])
    ylim(max_diff.*[-1 1])
    ax = gca; ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Difference in Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux (FBA - Algorithm)'], 'FontSize',titleSize)
    legend({'1 Model','2 Models, Model #1','2 Models, Model #2'}, 'FontSize',legendSize, 'Location','SouthEast', 'EdgeColor','w');
end

%% Constraint Violations
load('flux_color.mat')
tol = 1e-5;

if plotFlag == 3
    Sv_lim = max([max(abs(model_flux1{1}.warningFlag.Sv(:))), max(abs(model_flux2{1}.warningFlag.Sv(:))), max(abs(model_flux2{2}.warningFlag.Sv(:)))]);
    [mets_idx_1m,~] = find(abs(model_flux1{1}.warningFlag.Sv) > tol); mets_idx_1m = unique(mets_idx_1m);
    [mets_idx_2m1,~] = find(abs(model_flux2{1}.warningFlag.Sv) > tol); mets_idx_2m1 = unique(mets_idx_2m1);
    [mets_idx_2m2,~] = find(abs(model_flux2{2}.warningFlag.Sv) > tol); mets_idx_2m2 = unique(mets_idx_2m2);
    mets_idx = unique([mets_idx_1m; mets_idx_2m1; mets_idx_2m2]);
    
    [rxns_idx_1m,~] = find(model_flux1{1}.warningFlag.intl_bounds); rxns_idx_1m = unique(rxns_idx_1m);
    [rxns_idx_2m1,~] = find(model_flux2{1}.warningFlag.intl_bounds); rxns_idx_2m1 = unique(rxns_idx_2m1);
    [rxns_idx_2m2,~] = find(model_flux2{2}.warningFlag.intl_bounds); rxns_idx_2m2 = unique(rxns_idx_2m2);
    rxns_idx = unique([rxns_idx_1m; rxns_idx_2m1; rxns_idx_2m2]);
elseif plotFlag == 2
    Sv_lim = max([max(abs(model_flux2{1}.warningFlag.Sv(:))), max(abs(model_flux2{2}.warningFlag.Sv(:)))]);
    [mets_idx_2m1,~] = find(abs(model_flux2{1}.warningFlag.Sv) > tol); mets_idx_2m1 = unique(mets_idx_2m1);
    [mets_idx_2m2,~] = find(abs(model_flux2{2}.warningFlag.Sv) > tol); mets_idx_2m2 = unique(mets_idx_2m2);
    mets_idx = unique([mets_idx_2m1; mets_idx_2m2]);
    
    [rxns_idx_2m1,~] = find(model_flux2{1}.warningFlag.intl_bounds); rxns_idx_2m1 = unique(rxns_idx_2m1);
    [rxns_idx_2m2,~] = find(model_flux2{2}.warningFlag.intl_bounds); rxns_idx_2m2 = unique(rxns_idx_2m2);
    rxns_idx = unique([rxns_idx_2m1; rxns_idx_2m2]);
elseif plotFlag == 1
    Sv_lim = max(abs(model_flux1{1}.warningFlag.Sv(:)));
    [mets_idx_1m,~] = find(abs(model_flux1{1}.warningFlag.Sv) > tol);
    mets_idx = unique(mets_idx_1m);
    
    [rxns_idx_1m,~] = find(model_flux1{1}.warningFlag.intl_bounds);
    rxns_idx = unique(rxns_idx_1m);
end

if plotFlag == 3 || plotFlag == 1
    fig_handles(25) = figure(25);
    imagesc(model_flux1{1}.warningFlag.Sv(mets_idx,intlCon_idx1)); caxis(Sv_lim.*[-1 1]); colormap(flux_color)
    cbar = colorbar; cbar.FontSize = 0.75*axesLabelSize;
    cbar.Label.String = 'S\timesv'; cbar.Label.FontSize = xyLabelSize;
    ax = gca; ax.FontSize = 0.5*axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    ax.XTick = 1:numel(intlCon_idx1); ax.XTickLabel = model_flux1{1}.intl_con(intlCon_idx1);
    ax.YTick = 1:numel(mets_idx); ax.YTickLabel = model_flux1{1}.mets(mets_idx);
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Metabolites', 'FontSize',xyLabelSize)
    title([model_title ': 1 Model'], 'FontSize',titleSize)
end

if plotFlag == 3 || plotFlag == 2
    fig_handles(26) = figure(26);
    imagesc(model_flux2{1}.warningFlag.Sv(mets_idx,intlCon_idx2)); caxis(Sv_lim.*[-1 1]); colormap(flux_color)
    cbar = colorbar; cbar.FontSize = 0.75*axesLabelSize;
    cbar.Label.String = 'S\timesv'; cbar.Label.FontSize = xyLabelSize;
    ax = gca; ax.FontSize = 0.5*axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    ax.XTick = 1:numel(intlCon_idx2); ax.XTickLabel = model_flux2{1}.intl_con(intlCon_idx2);
    ax.YTick = 1:numel(mets_idx); ax.YTickLabel = model_flux2{1}.mets(mets_idx);
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Metabolites', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Model 1'], 'FontSize',titleSize)
    
    fig_handles(27) = figure(27);
    imagesc(model_flux2{2}.warningFlag.Sv(mets_idx,intlCon_idx2)); caxis(Sv_lim.*[-1 1]); colormap(flux_color)
    cbar = colorbar; cbar.FontSize = 0.75*axesLabelSize;
    cbar.Label.String = 'S\timesv'; cbar.Label.FontSize = xyLabelSize;
    ax = gca; ax.FontSize = 0.5*axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    ax.XTick = 1:numel(intlCon_idx2); ax.XTickLabel = model_flux2{2}.intl_con(intlCon_idx2);
    ax.YTick = 1:numel(mets_idx); ax.YTickLabel = model_flux2{2}.mets(mets_idx);
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Metabolites', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Model 2'], 'FontSize',titleSize)
end

if plotFlag == 3 || plotFlag == 1
    fig_handles(28) = figure(28);
    imagesc(model_flux1{1}.warningFlag.intl_bounds(rxns_idx,intlCon_idx1))
    caxis([0 1]); colormap(flip(gray(2))); cbar = colorbar; cbar.FontSize = 0.75*axesLabelSize;
    cbar.Ticks = [0.25 0.75]; cbar.TickLabels = [0,1];
    cbar.Label.String = 'Bounds Violation'; cbar.Label.FontSize = xyLabelSize;
    ax = gca; ax.FontSize = 0.5*axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    if size(model_flux1{1}.warningFlag.intl_bounds(rxns_idx,intlCon_idx1),1)~=0
        ax.XTick = 1:numel(intlCon_idx1); ax.XTickLabel = model_flux1{1}.intl_con(intlCon_idx1);
    else
        ax.XTick = [];
    end
    ax.YTick = 1:numel(rxns_idx); ax.YTickLabel = model_flux1{1}.rxns(model_flux1{1}.intl_idx(rxns_idx));
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Reactions', 'FontSize',xyLabelSize)
    title([model_title ': 1 Model'], 'FontSize',titleSize)
end

if plotFlag == 3 || plotFlag == 2
    fig_handles(29) = figure(29);
    imagesc(model_flux2{1}.warningFlag.intl_bounds(rxns_idx,intlCon_idx2))
    caxis([0 1]); colormap(flip(gray(2))); cbar = colorbar; cbar.FontSize = 0.75*axesLabelSize;
    cbar.Ticks = [0.25 0.75]; cbar.TickLabels = [0,1];
    cbar.Label.String = 'Bounds Violation'; cbar.Label.FontSize = xyLabelSize;
    ax = gca; ax.FontSize = 0.5*axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    if size(model_flux2{1}.warningFlag.intl_bounds(rxns_idx,intlCon_idx2),1)~=0
        ax.XTick = 1:numel(intlCon_idx2); ax.XTickLabel = model_flux2{1}.intl_con(intlCon_idx2);
    else
        ax.XTick = [];
    end
    ax.YTick = 1:numel(rxns_idx); ax.YTickLabel = model_flux2{1}.rxns(model_flux2{1}.intl_idx(rxns_idx));
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Reactions', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Model 1'], 'FontSize',titleSize)
    
    fig_handles(30) = figure(30);
    imagesc(model_flux2{2}.warningFlag.intl_bounds(rxns_idx,intlCon_idx2))
    caxis([0 1]); colormap(flip(gray(2))); cbar = colorbar; cbar.FontSize = 0.75*axesLabelSize;
    cbar.Ticks = [0.25 0.75]; cbar.TickLabels = [0,1];
    cbar.Label.String = 'Bounds Violation'; cbar.Label.FontSize = xyLabelSize;
    ax = gca; ax.FontSize = 0.5*axesLabelSize; ax.TickLength = tickLength.*[1 2];
    ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    if size(model_flux2{2}.warningFlag.intl_bounds(rxns_idx,intlCon_idx2),1)~=0
        ax.XTick = 1:numel(intlCon_idx2); ax.XTickLabel = model_flux2{2}.intl_con(intlCon_idx2);
    else
        ax.XTick = [];
    end
    ax.YTick = 1:numel(rxns_idx); ax.YTickLabel = model_flux2{2}.rxns(model_flux2{2}.intl_idx(rxns_idx));
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Reactions', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Model 2'], 'FontSize',titleSize)
end

%% Save Figures

if saveFlag
    if plotFlag == 3
        % Figure 1: Biomass Flux
        saveas(fig_handles(1),[plotFolder1 'biomassFlux'],'png');
        saveas(fig_handles(1),[plotFolder1 'biomassFlux'],'fig');
        saveas(fig_handles(1),[plotFolder2 'biomassFlux'],'png');
        saveas(fig_handles(1),[plotFolder2 'biomassFlux'],'fig');
        % Figure 2: Biomass Flux (Zoomed)
        saveas(fig_handles(2),[plotFolder1 'biomassFlux_zoomed'],'png');
        saveas(fig_handles(2),[plotFolder1 'biomassFlux_zoomed'],'fig');
        saveas(fig_handles(2),[plotFolder2 'biomassFlux_zoomed'],'png');
        saveas(fig_handles(2),[plotFolder2 'biomassFlux_zoomed'],'fig');
        % Figure 3: Intracellular Flux, 1 Model (Heat Map)
        saveas(fig_handles(3),[plotFolder1 'intlFlux_1m'],'png');
        saveas(fig_handles(3),[plotFolder1 'intlFlux_1m'],'fig');
        % Figure 4: Intracellular Flux, 2 Models, Model 1 (Heat Map)
        saveas(fig_handles(4),[plotFolder2 'intlFlux_2m1'],'png');
        saveas(fig_handles(4),[plotFolder2 'intlFlux_2m1'],'fig');
        % Figure 5: Intracellular Flux, 2 Models, Model 2 (Heat Map)
        saveas(fig_handles(5),[plotFolder2 'intlFlux_2m2'],'png');
        saveas(fig_handles(5),[plotFolder2 'intlFlux_2m2'],'fig');
        % Figure 6: Intracellular Binary Variables, 1 Model (Heat Map)
        saveas(fig_handles(6),[plotFolder1 'intlT_1m'],'png');
        saveas(fig_handles(6),[plotFolder1 'intlT_1m'],'fig');
        % Figure 7: Intracellular Binary Variables, 2 Models, Model 1 (Heat Map)
        saveas(fig_handles(7),[plotFolder2 'intlT_2m1'],'png');
        saveas(fig_handles(7),[plotFolder2 'intlT_2m1'],'fig');
        % Figure 8: Intracellular Binary Variables, 2 Models, Model 2 (Heat Map)
        saveas(fig_handles(8),[plotFolder2 'intlT_2m2'],'png');
        saveas(fig_handles(8),[plotFolder2 'intlT_2m2'],'fig');
        % Figure 9: 1 Model Clustered (Dendrogram)
        saveas(fig_handles(9),[plotFolder1 'dendrogram_1m'],'png');
        saveas(fig_handles(9),[plotFolder1 'dendrogram_1m'],'fig');
        % Figure 10: 2 Models, Model 1 Clustered (Dendrogram)
        saveas(fig_handles(10),[plotFolder2 'dendrogram_2m1'],'png');
        saveas(fig_handles(10),[plotFolder2 'dendrogram_2m1'],'fig');
        % Figure 11: 2 Models, Model 2 Clustered (Dendrogram)
        saveas(fig_handles(11),[plotFolder2 'dendrogram_2m2'],'png');
        saveas(fig_handles(11),[plotFolder2 'dendrogram_2m2'],'fig');
        % Figure 12: 1 Model & 2 Models Clustered (Dendrogram)
        saveas(fig_handles(12),[plotFolder1 'dendrogram'],'png');
        saveas(fig_handles(12),[plotFolder1 'dendrogram'],'fig');
        saveas(fig_handles(12),[plotFolder2 'dendrogram'],'png');
        saveas(fig_handles(12),[plotFolder2 'dendrogram'],'fig');
        % Figure 13: Intracellular Flux t-SNE
        saveas(fig_handles(13),[plotFolder1 'tSNE'],'png');
        saveas(fig_handles(13),[plotFolder1 'tSNE'],'fig');
        saveas(fig_handles(13),[plotFolder2 'tSNE'],'png');
        saveas(fig_handles(13),[plotFolder2 'tSNE'],'fig');
        % Figure 14: Intracellular Flux Correlation, 2 Models (Heat Map)
        saveas(fig_handles(14),[plotFolder2 'intlFluxCorr'],'png');
        saveas(fig_handles(14),[plotFolder2 'intlFluxCorr'],'fig');
        % Figure 15: Exchange Flux, 1 Model (Heat Map)
        saveas(fig_handles(15),[plotFolder1 'exchFlux_1m'],'png');
        saveas(fig_handles(15),[plotFolder1 'exchFlux_1m'],'fig');
        % Figure 16: Exchange Flux, 2 Models, Model 1 (Heat Map)
        saveas(fig_handles(16),[plotFolder2 'exchFlux_2m1'],'png');
        saveas(fig_handles(16),[plotFolder2 'exchFlux_2m1'],'fig');
        % Figure 17: Exchange Flux, 2 Models, Model 2 (Heat Map)
        saveas(fig_handles(17),[plotFolder2 'exchFlux_2m2'],'png');
        saveas(fig_handles(17),[plotFolder2 'exchFlux_2m2'],'fig');
        % Figure 18: Exchange Flux Correlation, 2 Models (Heat Map)
        saveas(fig_handles(18),[plotFolder2 'exchFluxCorr'],'png');
        saveas(fig_handles(18),[plotFolder2 'exchFluxCorr'],'fig');
        % Figure 19: Extracellular Metabolites that are Taken Up, Produced, or Exchanged (Heat Map)
        saveas(fig_handles(19),[plotFolder2 'exchMets_nonZeroFlux'],'png');
        saveas(fig_handles(19),[plotFolder2 'exchMets_nonZeroFlux'],'fig');
        % Figure 20: Metabolites that are Exchanged (Heat Map)
        saveas(fig_handles(20),[plotFolder2 'exchMets_exchanged'],'png');
        saveas(fig_handles(20),[plotFolder2 'exchMets_exchanged'],'fig');
        % Figure 21: Number of Times Metabolites are Exchanged
        saveas(fig_handles(21),[plotFolder2 'exchMets_exchangedNum'],'png');
        saveas(fig_handles(21),[plotFolder2 'exchMets_exchangedNum'],'fig');
        % Figure 22: Number of Times Metabolites are Exchanged by Intracellular at each Sparsity Constraint
        saveas(fig_handles(22),[plotFolder2 'exchMets_exchangedNum_intlCon'],'png');
        saveas(fig_handles(22),[plotFolder2 'exchMets_exchangedNum_intlCon'],'fig');
        % Figure 23: Biomass Flux Difference
        saveas(fig_handles(23),[plotFolder1 'biomassFluxDiff'],'png');
        saveas(fig_handles(23),[plotFolder1 'biomassFluxDiff'],'fig');
        saveas(fig_handles(23),[plotFolder2 'biomassFluxDiff'],'png');
        saveas(fig_handles(23),[plotFolder2 'biomassFluxDiff'],'fig');
        % Figure 24: Biomass Flux Difference (Zoomed)
        saveas(fig_handles(24),[plotFolder1 'biomassFluxDiff_zoomed'],'png');
        saveas(fig_handles(24),[plotFolder1 'biomassFluxDiff_zoomed'],'fig');
        saveas(fig_handles(24),[plotFolder2 'biomassFluxDiff_zoomed'],'png');
        saveas(fig_handles(24),[plotFolder2 'biomassFluxDiff_zoomed'],'fig');
        % Figure 25: Sv Violations, 1 Model
        saveas(fig_handles(25),[plotFolder1 'viol_Sv_1m'],'png');
        saveas(fig_handles(25),[plotFolder1 'viol_Sv_1m'],'fig');
        % Figure 26: Sv Violations, 2 Models, Model 1
        saveas(fig_handles(26),[plotFolder2 'viol_Sv_2m1'],'png');
        saveas(fig_handles(26),[plotFolder2 'viol_Sv_2m1'],'fig');
        % Figure 27: Sv Violations, 2 Models, Model 2
        saveas(fig_handles(27),[plotFolder2 'viol_Sv_2m2'],'png');
        saveas(fig_handles(27),[plotFolder2 'viol_Sv_2m2'],'fig');
        % Figure 28: Bounds Violations, 1 Model
        saveas(fig_handles(28),[plotFolder1 'viol_bounds_1m'],'png');
        saveas(fig_handles(28),[plotFolder1 'viol_bounds_1m'],'fig');
        % Figure 29: Bounds Violations, 2 Models, Model 1
        saveas(fig_handles(29),[plotFolder2 'viol_bounds_2m1'],'png');
        saveas(fig_handles(29),[plotFolder2 'viol_bounds_2m1'],'fig');
        % Figure 30: Bounds Violations, 2 Models, Model 2
        saveas(fig_handles(30),[plotFolder2 'viol_bounds_2m2'],'png');
        saveas(fig_handles(30),[plotFolder2 'viol_bounds_2m2'],'fig');
    elseif plotFlag == 2
        % Figure 1: Biomass Flux
        saveas(fig_handles(1),[plotFolder2 'biomassFlux'],'png');
        saveas(fig_handles(1),[plotFolder2 'biomassFlux'],'fig');
        % Figure 2: Biomass Flux (Zoomed)
        saveas(fig_handles(2),[plotFolder2 'biomassFlux_zoomed'],'png');
        saveas(fig_handles(2),[plotFolder2 'biomassFlux_zoomed'],'fig');
        % Figure 3: Intracellular Flux, 1 Model (Heat Map)
        % Figure 4: Intracellular Flux, 2 Models, Model 1 (Heat Map)
        saveas(fig_handles(4),[plotFolder2 'intlFlux_2m1'],'png');
        saveas(fig_handles(4),[plotFolder2 'intlFlux_2m1'],'fig');
        % Figure 5: Intracellular Flux, 2 Models, Model 2 (Heat Map)
        saveas(fig_handles(5),[plotFolder2 'intlFlux_2m2'],'png');
        saveas(fig_handles(5),[plotFolder2 'intlFlux_2m2'],'fig');
        % Figure 6: Intracellular Binary Variables, 1 Model (Heat Map)
        % Figure 7: Intracellular Binary Variables, 2 Models, Model 1 (Heat Map)
        saveas(fig_handles(7),[plotFolder2 'intlT_2m1'],'png');
        saveas(fig_handles(7),[plotFolder2 'intlT_2m1'],'fig');
        % Figure 8: Intracellular Binary Variables, 2 Models, Model 2 (Heat Map)
        saveas(fig_handles(8),[plotFolder2 'intlT_2m2'],'png');
        saveas(fig_handles(8),[plotFolder2 'intlT_2m2'],'fig');
        % Figure 9: 1 Model Clustered (Dendrogram)
        % Figure 10: 2 Models, Model 1 Clustered (Dendrogram)
        saveas(fig_handles(10),[plotFolder2 'dendrogram_2m1'],'png');
        saveas(fig_handles(10),[plotFolder2 'dendrogram_2m1'],'fig');
        % Figure 11: 2 Models, Model 2 Clustered (Dendrogram)
        saveas(fig_handles(11),[plotFolder2 'dendrogram_2m2'],'png');
        saveas(fig_handles(11),[plotFolder2 'dendrogram_2m2'],'fig');
        % Figure 12: 1 Model & 2 Models Clustered (Dendrogram)
        saveas(fig_handles(12),[plotFolder2 'dendrogram'],'png');
        saveas(fig_handles(12),[plotFolder2 'dendrogram'],'fig');
        % Figure 13: Intracellular Flux t-SNE
        saveas(fig_handles(13),[plotFolder2 'tSNE'],'png');
        saveas(fig_handles(13),[plotFolder2 'tSNE'],'fig');
        % Figure 14: Intracellular Flux Correlation, 2 Models (Heat Map)
        saveas(fig_handles(14),[plotFolder2 'intlFluxCorr'],'png');
        saveas(fig_handles(14),[plotFolder2 'intlFluxCorr'],'fig');
        % Figure 15: Exchange Flux, 1 Model (Heat Map)
        % Figure 16: Exchange Flux, 2 Models, Model 1 (Heat Map)
        saveas(fig_handles(16),[plotFolder2 'exchFlux_2m1'],'png');
        saveas(fig_handles(16),[plotFolder2 'exchFlux_2m1'],'fig');
        % Figure 17: Exchange Flux, 2 Models, Model 2 (Heat Map)
        saveas(fig_handles(17),[plotFolder2 'exchFlux_2m2'],'png');
        saveas(fig_handles(17),[plotFolder2 'exchFlux_2m2'],'fig');
        % Figure 18: Exchange Flux Correlation, 2 Models (Heat Map)
        saveas(fig_handles(18),[plotFolder2 'exchFluxCorr'],'png');
        saveas(fig_handles(18),[plotFolder2 'exchFluxCorr'],'fig');
        % Figure 19: Extracellular Metabolites that are Taken Up, Produced, or Exchanged (Heat Map)
        saveas(fig_handles(19),[plotFolder2 'exchMets_nonZeroFlux'],'png');
        saveas(fig_handles(19),[plotFolder2 'exchMets_nonZeroFlux'],'fig');
        % Figure 20: Metabolites that are Exchanged (Heat Map)
        saveas(fig_handles(20),[plotFolder2 'exchMets_exchanged'],'png');
        saveas(fig_handles(20),[plotFolder2 'exchMets_exchanged'],'fig');
        % Figure 21: Number of Times Metabolites are Exchanged
        saveas(fig_handles(21),[plotFolder2 'exchMets_exchangedNum'],'png');
        saveas(fig_handles(21),[plotFolder2 'exchMets_exchangedNum'],'fig');
        % Figure 22: Number of Times Metabolites are Exchanged by Intracellular at each Sparsity Constraint
        saveas(fig_handles(22),[plotFolder2 'exchMets_exchangedNum_intlCon'],'png');
        saveas(fig_handles(22),[plotFolder2 'exchMets_exchangedNum_intlCon'],'fig');
        % Figure 23: Biomass Flux Difference
        saveas(fig_handles(23),[plotFolder2 'biomassFluxDiff'],'png');
        saveas(fig_handles(23),[plotFolder2 'biomassFluxDiff'],'fig');
        % Figure 24: Biomass Flux Difference (Zoomed)
        saveas(fig_handles(24),[plotFolder2 'biomassFluxDiff_zoomed'],'png');
        saveas(fig_handles(24),[plotFolder2 'biomassFluxDiff_zoomed'],'fig');
        % Figure 25: Sv Violations, 1 Model
        % Figure 26: Sv Violations, 2 Models, Model 1
        saveas(fig_handles(26),[plotFolder2 'viol_Sv_2m1'],'png');
        saveas(fig_handles(26),[plotFolder2 'viol_Sv_2m1'],'fig');
        % Figure 27: Sv Violations, 2 Models, Model 2
        saveas(fig_handles(27),[plotFolder2 'viol_Sv_2m2'],'png');
        saveas(fig_handles(27),[plotFolder2 'viol_Sv_2m2'],'fig');
        % Figure 28: Bounds Violations, 1 Model
        % Figure 29: Bounds Violations, 2 Models, Model 1
        saveas(fig_handles(29),[plotFolder2 'viol_bounds_2m1'],'png');
        saveas(fig_handles(29),[plotFolder2 'viol_bounds_2m1'],'fig');
        % Figure 30: Bounds Violations, 2 Models, Model 2
        saveas(fig_handles(30),[plotFolder2 'viol_bounds_2m2'],'png');
        saveas(fig_handles(30),[plotFolder2 'viol_bounds_2m2'],'fig');
    elseif plotFlag == 1
        % Figure 1: Biomass Flux
        saveas(fig_handles(1),[plotFolder1 'biomassFlux'],'png');
        saveas(fig_handles(1),[plotFolder1 'biomassFlux'],'fig');
        % Figure 2: Biomass Flux (Zoomed)
        saveas(fig_handles(2),[plotFolder1 'biomassFlux_zoomed'],'png');
        saveas(fig_handles(2),[plotFolder1 'biomassFlux_zoomed'],'fig');
        % Figure 3: Intracellular Flux, 1 Model (Heat Map)
        saveas(fig_handles(3),[plotFolder1 'intlFlux_1m'],'png');
        saveas(fig_handles(3),[plotFolder1 'intlFlux_1m'],'fig');
        % Figure 4: Intracellular Flux, 2 Models, Model 1 (Heat Map)
        % Figure 5: Intracellular Flux, 2 Models, Model 2 (Heat Map)
        % Figure 6: Intracellular Binary Variables, 1 Model (Heat Map)
        saveas(fig_handles(6),[plotFolder1 'intlT_1m'],'png');
        saveas(fig_handles(6),[plotFolder1 'intlT_1m'],'fig');
        % Figure 7: Intracellular Binary Variables, 2 Models, Model 1 (Heat Map)
        % Figure 8: Intracellular Binary Variables, 2 Models, Model 2 (Heat Map)
        % Figure 9: 1 Model Clustered (Dendrogram)
        saveas(fig_handles(9),[plotFolder1 'dendrogram_1m'],'png');
        saveas(fig_handles(9),[plotFolder1 'dendrogram_1m'],'fig');
        % Figure 10: 2 Models, Model 1 Clustered (Dendrogram)
        % Figure 11: 2 Models, Model 2 Clustered (Dendrogram)
        % Figure 12: 1 Model & 2 Models Clustered (Dendrogram)
        saveas(fig_handles(12),[plotFolder1 'dendrogram'],'png');
        saveas(fig_handles(12),[plotFolder1 'dendrogram'],'fig');
        % Figure 13: Intracellular Flux t-SNE
        % Figure 14: Intracellular Flux Correlation, 2 Models (Heat Map)
        % Figure 15: Exchange Flux, 1 Model (Heat Map)
        saveas(fig_handles(15),[plotFolder1 'exchFlux_1m'],'png');
        saveas(fig_handles(15),[plotFolder1 'exchFlux_1m'],'fig');
        % Figure 16: Exchange Flux, 2 Models, Model 1 (Heat Map)
        % Figure 17: Exchange Flux, 2 Models, Model 2 (Heat Map)
        % Figure 18: Exchange Flux Correlation, 2 Models (Heat Map)
        % Figure 19: Extracellular Metabolites that are Taken Up, Produced, or Exchanged (Heat Map)
        % Figure 20: Metabolites that are Exchanged (Heat Map)
        % Figure 21: Number of Times Metabolites are Exchanged
        % Figure 22: Number of Times Metabolites are Exchanged by Intracellular at each Sparsity Constraint
        % Figure 23: Biomass Flux Difference
        saveas(fig_handles(23),[plotFolder1 'biomassFluxDiff'],'png');
        saveas(fig_handles(23),[plotFolder1 'biomassFluxDiff'],'fig');
        % Figure 24: Biomass Flux Difference (Zoomed)
        saveas(fig_handles(24),[plotFolder1 'biomassFluxDiff_zoomed'],'png');
        saveas(fig_handles(24),[plotFolder1 'biomassFluxDiff_zoomed'],'fig');
        % Figure 25: Sv Violations, 1 Model
        saveas(fig_handles(25),[plotFolder1 'viol_Sv_1m'],'png');
        saveas(fig_handles(25),[plotFolder1 'viol_Sv_1m'],'fig');
        % Figure 26: Sv Violations, 2 Models, Model 1
        % Figure 27: Sv Violations, 2 Models, Model 2
        % Figure 28: Bounds Violations, 1 Model
        saveas(fig_handles(28),[plotFolder1 'viol_bounds_1m'],'png');
        saveas(fig_handles(28),[plotFolder1 'viol_bounds_1m'],'fig');
        % Figure 29: Bounds Violations, 2 Models, Model 1
        % Figure 30: Bounds Violations, 2 Models, Model 2
    end
end

end

% Bi-Cluster Flux
function [P1,P2] = biclusterFlux(flux_data,method,metric)
if ~exist('method','var')
    method = 'single';
end
if ~exist('metric','var')
    metric = 'euclidean';
end

% Cluster by Columns
Z1 = linkage(double(flux_data),method,metric); % get tree of hierarchical clusters (cluster closest distances)
% Cluster by Rows
Z2 = linkage(double(flux_data)',method,metric); % get tree of hierarchical clusters (cluster closest distances)
% Dendrogram
figure; [~, ~, P1] = dendrogram(Z1,0); close gcf
figure; [~, ~, P2] = dendrogram(Z2,0); close gcf
end