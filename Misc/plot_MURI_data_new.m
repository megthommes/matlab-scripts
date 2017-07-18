function plot_MURI_data_new(model1, model2, plotFolder, model_title, medium_rxns)
%plot_MURI_data_new Create all the plots for the MURI data
%   plot_MURI_data_new(model1, model2, plotFolder, model_title)
%
%INPUTS
% model*: cell structure that contains the following fields:
%   sparse_con: sparsity constraints, vector
%   mets: metabolite names, vector
%   rxns: reaction names, vector
%   biomass: biomass flux, vector
%   flux: reaction fluxes, matrix [rxns x sparse_con]
%   int: reaction binary variables t, matrix [rxns x sparse_con]
%   exch_idx: index of exchange ("uptake") reactions
%   trspt_idx: index of transport ("exchange") reactions
%   intl_idx: index of intracellular ("internal") reactions
% plotFolder: cell array path to save figures. First index is for model1,
%   second index is for model2
%
% model_title: name of model
% medium_rxns: name of reactions in the medium
%
% If you are not plotting both models (model1 and model2), then enter "''"
% for the model you are not plotting.
%   Example: plot_MURI_data_new(model1, '', plotFolder, model_title)
%
% Meghan Thommes 07/13/2017

axesLabelSize = 14;
xyLabelSize = 18;
titleSize = 22;
legendSize = 8;
lineWidth = 3;
arrowHead = 2;
markerSize = 10;

if ~isempty(model1) && ~isempty(model2)
    plotFolder1 = plotFolder{1};
    plotFolder2 = plotFolder{2};
    sparseCon_idx1 = find(model1{1}.biomass~=0);
    sparseCon_idx2 = find(model2{1}.biomass~=0);
    plotFlag = 3;
elseif isempty(model1)
    plotFolder2 = plotFolder{2};
    sparseCon_idx2 = find(model2{1}.biomass~=0);
    plotFlag = 2;
elseif isempty(model2)
    plotFolder1 = plotFolder{1};
    sparseCon_idx1 = find(model1{1}.biomass~=0);
    plotFlag = 1;
end

if exist('medium_rxns','var')
    medium_rxns = strrep(strrep(medium_rxns,'EX_',''),'_e','');
end

color = parula(3);
%% Biomass

if plotFlag == 3 % model1 and model2
    figure(1);
    plot(model1{1}.sparse_con,model1{1}.biomass,'ro-', 'LineWidth',lineWidth); hold on
    plot(model2{1}.sparse_con,model2{1}.biomass+model2{2}.biomass,'o--', 'Color',[256 150 150]./256, 'LineWidth',lineWidth);
    plot(model2{1}.sparse_con,model2{1}.biomass,'gd-', 'LineWidth',lineWidth);
    plot(model2{2}.sparse_con,model2{2}.biomass,'bs--', 'LineWidth',lineWidth); hold off
    xlim([min([model1{1}.sparse_con, model2{1}.sparse_con]) max([model1{1}.sparse_con, model2{1}.sparse_con])])
    set(gca, 'FontSize',axesLabelSize); grid on
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux'], 'FontSize',titleSize)
    l_h = legend({'Biomass of 1 Model','Total Biomass of 2 Models','Biomass of Model #1','Biomass of Model #2'}, 'Location','Best'); set(l_h, 'EdgeColor','w');
    saveas(gcf,[plotFolder1 'biomassFlux'],'png'); saveas(gcf,[plotFolder2 'biomassFlux'],'png');
    saveas(gcf,[plotFolder1 'biomassFlux'],'fig'); saveas(gcf,[plotFolder2 'biomassFlux'],'fig');
    
    model1_max = model1{1}.sparse_con(sparseCon_idx1(find(abs(diff(model1{1}.biomass(sparseCon_idx1))) <= 1E-9,1))); if isempty(model1_max); model1_max = NaN; end
    model2_max = model2{1}.sparse_con(sparseCon_idx2(find(abs(diff(model2{1}.biomass(sparseCon_idx2))) <= 1E-9,1))); if isempty(model2_max); model2_max = NaN; end
    model1_min = model1{1}.sparse_con(find(model1{1}.biomass == 0, 1, 'last')); if isempty(model1_min); model1_min = NaN; end
    model2_min = model2{1}.sparse_con(find(model2{1}.biomass == 0, 1, 'last')); if isempty(model2_min); model2_min = NaN; end
    figure(2);
    plot(model1{1}.sparse_con,model1{1}.biomass,'ro-', 'LineWidth',lineWidth); hold on
    plot(model2{1}.sparse_con,model2{1}.biomass+model2{2}.biomass,'o--', 'Color',[256 150 150]./256, 'LineWidth',lineWidth);
    plot(model2{1}.sparse_con,model2{1}.biomass,'gd-', 'LineWidth',lineWidth);
    plot(model2{2}.sparse_con,model2{2}.biomass,'bs--', 'LineWidth',lineWidth); hold off
    set(gca, 'FontSize',axesLabelSize); grid on
    xlim([nanmin(model1_min,model2_min)-1 nanmax(model1_max,model2_max)+1])
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux'], 'FontSize',titleSize)
    l_h = legend({'Biomass of 1 Model','Total Biomass of 2 Models','Biomass of Model #1','Biomass of Model #2'}, 'Location','Best'); set(l_h, 'EdgeColor','w');
    saveas(gcf,[plotFolder1 'biomassFlux_zoomed'],'png'); saveas(gcf,[plotFolder2 'biomassFlux_zoomed'],'png');
    saveas(gcf,[plotFolder1 'biomassFlux_zoomed'],'fig'); saveas(gcf,[plotFolder2 'biomassFlux_zoomed'],'fig');
    
elseif plotFlag == 2 % model2 only
    figure(1);
    plot(model2{1}.sparse_con,model2{1}.biomass+model2{2}.biomass,'o--', 'Color',[256 150 150]./256, 'LineWidth',lineWidth); hold on
    plot(model2{1}.sparse_con,model2{1}.biomass,'gd-', 'LineWidth',lineWidth);
    plot(model2{2}.sparse_con,model2{2}.biomass,'bs--', 'LineWidth',lineWidth); hold off
    xlim([min(model2{1}.sparse_con) max(model2{1}.sparse_con)])
    set(gca, 'FontSize',axesLabelSize); grid on
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux'], 'FontSize',titleSize)
    l_h = legend({'Total Biomass of 2 Models','Biomass of Model #1','Biomass of Model #2'}, 'Location','Best'); set(l_h, 'EdgeColor','w');
    saveas(gcf,[plotFolder2 'biomassFlux'],'png');
    saveas(gcf,[plotFolder2 'biomassFlux'],'fig');
    
    model2_max = model2{1}.sparse_con(sparseCon_idx2(find(abs(diff(model2{1}.biomass(sparseCon_idx2))) <= 1E-9,1))); if isempty(model2_max); model2_max = NaN; end
    model2_min = model2{1}.sparse_con(find(model2{1}.biomass == 0, 1, 'last')); if isempty(model2_min); model2_min = NaN; end
    figure(2);
    plot(model2{1}.sparse_con,model2{1}.biomass+model2{2}.biomass,'o--', 'Color',[256 150 150]./256, 'LineWidth',lineWidth); hold on
    plot(model2{1}.sparse_con,model2{1}.biomass,'gd-', 'LineWidth',lineWidth);
    plot(model2{2}.sparse_con,model2{2}.biomass,'bs--', 'LineWidth',lineWidth); hold off
    set(gca, 'FontSize',axesLabelSize); grid on
    xlim([nanmin(model2_min)-1 nanmax(model2_max)+1])
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux'], 'FontSize',titleSize)
    l_h = legend({'Total Biomass of 2 Models','Biomass of Model #1','Biomass of Model #2'}, 'Location','Best'); set(l_h, 'EdgeColor','w');
    saveas(gcf,[plotFolder2 'biomassFlux_zoomed'],'png');
    saveas(gcf,[plotFolder2 'biomassFlux_zoomed'],'fig');
    
elseif plotFlag == 1 % model1 only
    figure(1);
    plot(model1{1}.sparse_con,model1{1}.biomass,'ro-', 'LineWidth',lineWidth);
    xlim([min(model1{1}.sparse_con) max(model1{1}.sparse_con)])
    set(gca, 'FontSize',axesLabelSize); grid on
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux'], 'FontSize',titleSize)
    l_h = legend({'Biomass of 1 Model'}, 'Location','Best'); set(l_h, 'EdgeColor','w');
    saveas(gcf,[plotFolder1 'biomassFlux'],'png');
    saveas(gcf,[plotFolder1 'biomassFlux'],'fig');
    
    model1_max = model1{1}.sparse_con(sparseCon_idx1(find(abs(diff(model1{1}.biomass(sparseCon_idx1))) <= 1E-9,1))); if isempty(model1_max); model1_max = NaN; end
    model1_min = model1{1}.sparse_con(find(model1{1}.biomass == 0, 1, 'last')); if isempty(model1_min); model1_min = NaN; end
    figure(2);
    plot(model1{1}.sparse_con,model1{1}.biomass,'ro-', 'LineWidth',lineWidth);
    set(gca, 'FontSize',axesLabelSize); grid on
    xlim([nanmin(model1_min)-1 nanmax(model1_max)+1])
    xlabel('Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel('Biomass Flux (hr^{-1})', 'FontSize',xyLabelSize)
    title([model_title ': Biomass Flux'], 'FontSize',titleSize)
    l_h = legend({'Biomass of 1 Model'}, 'Location','Best'); set(l_h, 'EdgeColor','w');
    saveas(gcf,[plotFolder1 'biomassFlux_zoomed'],'png');
    saveas(gcf,[plotFolder1 'biomassFlux_zoomed'],'fig');
    
end
%% Exchange Type

if plotFlag == 3 || plotFlag == 2 % model1 and model2 or model2 only
    % Metabolites that are Taken Up, Produced, or Exchanged
    fig = figure(3); ax = gca;
    [idx_nonZero,~] = find(model2{1}.fluxType(:,sparseCon_idx2) ~= 0); idx_nonZero = unique(idx_nonZero);
    [fig,ax,cbar] = plotExchangeType(model2{1}.fluxType(idx_nonZero,sparseCon_idx2)',fig,ax,strrep(strrep(model2{1}.rxns(model2{1}.exch_idx(idx_nonZero)),'EX_',''),'_e',''));
    ax.XTick = 1:numel(sparseCon_idx2); ax.XTickLabel = model2{1}.sparse_con(sparseCon_idx2);
    ax.FontSize = 0.75*legendSize; ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    cbar.FontSize = axesLabelSize; cbar.Label.FontSize = xyLabelSize;
    xlabel(ax, 'Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel(ax, 'Extracellular Metabolites', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Metabolite Exchange'], 'FontSize',titleSize)
    box on; grid on
    saveas(gcf,[plotFolder2 'exchType_nonZero'],'png');
    saveas(gcf,[plotFolder2 'exchType_nonZero'],'fig');
    
    % Metabolites that are Exchanged
    fig = figure(4); ax = gca;
    [idx_7or8,~] = find(model2{1}.fluxType(:,sparseCon_idx2(1:end-1)) == 7 | model2{1}.fluxType(:,sparseCon_idx2(1:end-1)) == 8); idx_7or8 = unique(idx_7or8);
    [fig,ax,cbar] = plotExchangeType(model2{1}.fluxType(idx_7or8,sparseCon_idx2)',fig,ax,strrep(strrep(model2{1}.rxns(model2{1}.exch_idx(idx_7or8)),'EX_',''),'_e',''));
    ax.XTick = 1:numel(sparseCon_idx2); ax.XTickLabel = model2{1}.sparse_con(sparseCon_idx2);
    ax.FontSize = legendSize; ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    cbar.FontSize = axesLabelSize; cbar.Label.FontSize = xyLabelSize;
    xlabel(ax, 'Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel(ax, 'Extracellular Metabolites', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Metabolite Exchange'], 'FontSize',titleSize)
    box on; grid on
    saveas(gcf,[plotFolder2 'exchType_7or8'],'png');
    saveas(gcf,[plotFolder2 'exchType_7or8'],'fig');
    
    % Metabolites that are Exchanged and Not in the Medium
    fig = figure(5); ax = gca;
    [~,rm_idx,~] = intersect(strrep(strrep(model2{1}.rxns(model2{1}.exch_idx(idx_7or8)),'EX_',''),'_e',''),medium_rxns);
    idx_7or8_notInMedium = idx_7or8; idx_7or8_notInMedium(rm_idx) = [];
    [fig,ax,cbar] = plotExchangeType(model2{1}.fluxType(idx_7or8_notInMedium,sparseCon_idx2)',fig,ax,strrep(strrep(model2{1}.rxns(model2{1}.exch_idx(idx_7or8_notInMedium)),'EX_',''),'_e',''));
    ax.XTick = 1:numel(sparseCon_idx2); ax.XTickLabel = model2{1}.sparse_con(sparseCon_idx2);
    ax.FontSize = legendSize; ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 45;
    cbar.FontSize = axesLabelSize; cbar.Label.FontSize = xyLabelSize;
    xlabel(ax, 'Constraint Bound for the Number of Intracellular Reactions', 'FontSize',xyLabelSize)
    ylabel(ax, 'Extracellular Metabolites', 'FontSize',xyLabelSize)
    title([model_title ': 2 Models, Metabolite Exchange'], 'FontSize',titleSize)
    box on; grid on
    saveas(gcf,[plotFolder2 'exchType_7or8_notInMedium'],'png');
    saveas(gcf,[plotFolder2 'exchType_7or8_notInMedium'],'fig');
end

%% Crossfeedogram

if plotFlag == 3 || plotFlag == 2 % model1 and model2 or model2 only
    % Metabolites that are Exchanged
    fig = figure(6);
    for met_idx = 1:numel(idx_7or8)
        model1_flux = fliplr(model2{1}.flux(model2{1}.exch_idx(idx_7or8(met_idx)),sparseCon_idx2));
        model2_flux = fliplr(model2{2}.flux(model2{1}.exch_idx(idx_7or8(met_idx)),sparseCon_idx2));
        
        ax = subplot(round(sqrt(numel(idx_7or8))),ceil(numel(idx_7or8)/round(sqrt(numel(idx_7or8)))),met_idx);
        plotCocultureExchange(model1_flux,model2_flux,1,fig,ax,markerSize,...
            arrowHead,0.25*lineWidth,[1,0,0]);
        box on; grid on; hold off
        
        set(gca, 'FontSize',0.5*axesLabelSize)
        title(strrep(strrep(model2{1}.rxns(model2{1}.exch_idx(idx_7or8(met_idx))),'EX_',''),'_e',''), 'FontSize',0.5*xyLabelSize, 'Interpreter','none'); xlabel(ax,''); ylabel(ax,'')
    end
    [~,h] = suplabel('Model 1 Exchange Flux (mmol/gCDW/hr)'); set(h, 'FontSize',0.8*titleSize, 'Position',[0.5,0,0]);
    [~,h] = suplabel('Model 2 Exchange Flux (mmol/gCDW/hr)','y'); set(h, 'FontSize',0.8*titleSize, 'Position',[0,0.5,0])
    [~,h] = suplabel([model_title ': 2 Models (Co), Crossfeedogram'],'t'); set(h, 'FontSize',titleSize);
    saveas(fig,[plotFolder2 'crossfeedogram'],'png');
    saveas(fig,[plotFolder2 'crossfeedogram'],'fig');
    
    % Metabolites that are Exchanged and Not in the Medium
    fig = figure(7);
    for met_idx = 1:numel(idx_7or8_notInMedium)
        model1_flux = fliplr(model2{1}.flux(model2{1}.exch_idx(idx_7or8_notInMedium(met_idx)),sparseCon_idx2));
        model2_flux = fliplr(model2{2}.flux(model2{1}.exch_idx(idx_7or8_notInMedium(met_idx)),sparseCon_idx2));
        
        ax = subplot(round(sqrt(numel(idx_7or8_notInMedium))),ceil(numel(idx_7or8_notInMedium)/round(sqrt(numel(idx_7or8_notInMedium)))),met_idx);
        plotCocultureExchange(model1_flux,model2_flux,1,fig,ax,markerSize,...
            arrowHead,0.25*lineWidth,[1,0,0]);
        box on; grid on; hold off
        
        set(gca, 'FontSize',0.5*axesLabelSize)
        title(strrep(strrep(model2{1}.rxns(model2{1}.exch_idx(idx_7or8_notInMedium(met_idx))),'EX_',''),'_e',''), 'FontSize',0.5*xyLabelSize, 'Interpreter','none'); xlabel(ax,''); ylabel(ax,'')
    end
    [~,h] = suplabel('Model 1 Exchange Flux (mmol/gCDW/hr)'); set(h, 'FontSize',0.8*titleSize, 'Position',[0.5,0,0]);
    [~,h] = suplabel('Model 2 Exchange Flux (mmol/gCDW/hr)','y'); set(h, 'FontSize',0.8*titleSize, 'Position',[0,0.5,0])
    [~,h] = suplabel([model_title ': 2 Models (Co), Crossfeedogram'],'t'); set(h, 'FontSize',titleSize);
    saveas(fig,[plotFolder2 'crossfeedogram_notInMedium'],'png');
    saveas(fig,[plotFolder2 'crossfeedogram_notInMedium'],'fig');
end
%% Exchanged Metabolites

if plotFlag == 3 || plotFlag == 2 % model1 and model2 or model2 only
    headers = {'1->2','2->1'};
    col_num = 1;
    if exist([plotFolder2 'exchanged_metabolites.xlsx'],'file') ~= 0; delete([plotFolder2 'exchanged_metabolites.xlsx']); end
    model2{1}.exchMets_12 = cell(size(model2{1}.sparse_con));
    model2{1}.exchMets_21 = cell(size(model2{1}.sparse_con));
    for ii = 1:numel(sparseCon_idx2)
        [idx_12,~] = find(model2{1}.fluxType(:,sparseCon_idx2(ii)) == 8); idx_12 = unique(idx_12);
        exchMets_12 = strrep(strrep(model2{1}.rxns(model2{1}.exch_idx(idx_12)),'EX_',''),'_e','');
        model2{1}.exchMets_12{sparseCon_idx2(ii)} = exchMets_12;
        
        [idx_21,~] = find(model2{1}.fluxType(:,sparseCon_idx2(ii)) == 7); idx_21 = unique(idx_21);
        exchMets_21 = strrep(strrep(model2{1}.rxns(model2{1}.exch_idx(idx_21)),'EX_',''),'_e','');
        model2{1}.exchMets_21{sparseCon_idx2(ii)} = exchMets_21;
        
        exch_12 = [{int2str(model2{1}.sparse_con(sparseCon_idx2(ii)))}; headers{1}; exchMets_12];
        exch_21 = [{int2str(model2{1}.sparse_con(sparseCon_idx2(ii)))}; headers{2}; exchMets_21];
        writetable(table(exch_12),[plotFolder2 'exchanged_metabolites.xlsx'], 'WriteVariableNames',false, 'Range',[num2letter(col_num) '1']);
        writetable(table(exch_21),[plotFolder2 'exchanged_metabolites.xlsx'], 'WriteVariableNames',false, 'Range',[num2letter(col_num+1) '1'])
        col_num = col_num+2;
    end
    model2{2}.exchMets_12 = model2{1}.exchMets_12;
    model2{2}.exchMets_21 = model2{1}.exchMets_21;
    
    exchMets = union(cat(1, model2{1}.exchMets_12{:}), cat(1, model2{1}.exchMets_21{:}));
    model2{1}.exchMets = exchMets;
    model2{2}.exchMets = exchMets;
    
    num_12 = arrayfun(@(y) sum(arrayfun(@(x) sum(ismember(model2{1}.exchMets_12{x},exchMets{y})),1:numel(model2{1}.sparse_con))), 1:numel(exchMets));
    num_21 = arrayfun(@(y) sum(arrayfun(@(x) sum(ismember(model2{1}.exchMets_21{x},exchMets{y})),1:numel(model2{1}.sparse_con))), 1:numel(exchMets));
    num_total = num_12 + num_21;
    
    figure(8)
    bar([num_12', num_21', num_total'],1)
    xlim([0,numel(exchMets)]+0.5); ylim([0 max(num_total)+0.5]); box on; grid on
    set(gca, 'XTick',1:numel(exchMets), 'XTickLabel',exchMets, ...
        'XTickLabelRotation',45,'FontSize',axesLabelSize)
    xlabel('Extracellular Metabolites', 'FontSize',xyLabelSize)
    ylabel('Number of Times Exchanged', 'FontSize',xyLabelSize)
    title([model_title ': Exchanged Metabolites'], 'FontSize',titleSize)
    l_h = legend([headers, 'total'], 'FontSize',axesLabelSize, 'Location','NorthEast'); set(l_h, 'EdgeColor','w');
    saveas(gcf,[plotFolder2 'exchMets'],'png');
    saveas(gcf,[plotFolder2 'exchMets'],'fig');
    
    figure(9);
    if max(num_total) > 20
        edges = 0:ceil(max(num_total)/15):max(num_total);
        if max(num_total) > edges(end)
            edges(end+1) = edges(end) + ceil(max(num_total)/15);
        end
        edges = edges - ceil(max(num_total)/15)/2;
    else
        edges = 0:max(num_total)+1;
        edges = edges - 0.5;
    end
    h12 = histogram(num_12,edges, 'FaceColor', color(1,:), 'FaceAlpha',0.75); hold on
    h21 = histogram(num_21,edges, 'FaceColor', color(2,:), 'FaceAlpha',0.75); hold off
    xlim([edges(1) edges(end)]); ylim([0 max([h12.Values, h21.Values])+0.5]); box on; grid on
    set(gca, 'YTick',unique(ceil(get(gca, 'YTick'))), 'XTick',unique(ceil(get(gca, 'XTick'))), 'FontSize',axesLabelSize);
    xlabel('Number of Times Exchanged', 'FontSize',xyLabelSize)
    ylabel('Number of Extracellular Metabolites', 'FontSize',xyLabelSize)
    title([model_title ': Exchanged Metabolites'], 'FontSize',titleSize)
    l_h = legend(headers, 'FontSize',axesLabelSize, 'Location','NorthEast'); set(l_h, 'EdgeColor','w');
    saveas(gcf,[plotFolder2 'exchMets_hist'],'png');
    saveas(gcf,[plotFolder2 'exchMets_hist'],'fig');
end
%% Exchanged Metabolites that are not in the Medium

if plotFlag == 3 || plotFlag == 2 % model1 and model2 or model2 only
    col_num = 1;
    if exist([plotFolder2 'exchanged_metabolites_notInMedium.xlsx'],'file') ~= 0; delete([plotFolder2 'exchanged_metabolites_notInMedium.xlsx']); end
    model2{1}.exchMets_notInMedium_12 = cell(size(model2{1}.sparse_con));
    model2{1}.exchMets_notInMedium_21 = cell(size(model2{1}.sparse_con));
    for ii = 1:numel(sparseCon_idx2)-1
        exchMets_12 = model2{1}.exchMets_12{sparseCon_idx2(ii)};
        [~,rm_idx,~] = intersect(exchMets_12,medium_rxns); exchMets_12(rm_idx) = []; clear rm_idx
        model2{1}.exchMets_notInMedium_12{sparseCon_idx2(ii)} = exchMets_12;
        
        exchMets_21 = model2{1}.exchMets_21{sparseCon_idx2(ii)};
        [~,rm_idx,~] = intersect(exchMets_21,medium_rxns); exchMets_21(rm_idx) = []; clear rm_idx
        model2{1}.exchMets_notInMedium_21{sparseCon_idx2(ii)} = exchMets_21;
        
        exch_12 = [{int2str(model2{1}.sparse_con(sparseCon_idx2(ii)))}; headers{1}; exchMets_12];
        exch_21 = [{int2str(model2{1}.sparse_con(sparseCon_idx2(ii)))}; headers{2}; exchMets_21];
        writetable(table(exch_12),[plotFolder2 'exchanged_metabolites_notInMedium.xlsx'], 'WriteVariableNames',false, 'Range',[num2letter(col_num) '1']);
        writetable(table(exch_21),[plotFolder2 'exchanged_metabolites_notInMedium.xlsx'], 'WriteVariableNames',false, 'Range',[num2letter(col_num+1) '1'])
        col_num = col_num+2; clear idx* exchMets*
    end
    model2{2}.exchMets_notInMedium_12 = model2{1}.exchMets_notInMedium_12;
    model2{2}.exchMets_notInMedium_21 = model2{1}.exchMets_notInMedium_21;
    
    exchMets = union(cat(1, model2{1}.exchMets_notInMedium_12{:}), cat(1, model2{1}.exchMets_notInMedium_21{:}));
    model2{1}.exchMets_notInMedium = exchMets;
    model2{2}.exchMets_notInMedium = exchMets;
    
    num_12 = arrayfun(@(y) sum(arrayfun(@(x) sum(ismember(model2{1}.exchMets_notInMedium_12{x},exchMets{y})),1:numel(model2{1}.sparse_con))), 1:numel(exchMets));
    num_21 = arrayfun(@(y) sum(arrayfun(@(x) sum(ismember(model2{1}.exchMets_notInMedium_21{x},exchMets{y})),1:numel(model2{1}.sparse_con))), 1:numel(exchMets));
    num_total = num_12 + num_21;
    
    figure(10)
    bar([num_12', num_21', num_total'],1)
    xlim([0,numel(exchMets)]+0.5); ylim([0 max(num_total)+0.5]); box on; grid on
    set(gca, 'XTick',1:numel(exchMets), 'XTickLabel',exchMets, ...
        'XTickLabelRotation',45,'FontSize',axesLabelSize)
    xlabel('Extracellular Metabolites', 'FontSize',xyLabelSize)
    ylabel('Number of Times Exchanged', 'FontSize',xyLabelSize)
    title([model_title ': Exchanged Metabolites'], 'FontSize',titleSize)
    l_h = legend([headers, 'total'], 'FontSize',axesLabelSize, 'Location','NorthEast'); set(l_h, 'EdgeColor','w');
    saveas(gcf,[plotFolder2 'exchMets_notInMedium'],'png');
    saveas(gcf,[plotFolder2 'exchMets_notInMedium'],'fig');
    
    figure(11);
    if max(num_total) > 20
        edges = 0:ceil(max(num_total)/15):max(num_total);
        if max(num_total) > edges(end)
            edges(end+1) = edges(end) + ceil(max(num_total)/15);
        end
        edges = edges - ceil(max(num_total)/15)/2;
    else
        edges = 0:max(num_total)+1;
        edges = edges - 0.5;
    end
    h12 = histogram(num_12,edges, 'FaceColor', color(1,:), 'FaceAlpha',0.75); hold on
    h21 = histogram(num_21,edges, 'FaceColor', color(2,:), 'FaceAlpha',0.75); hold off
    xlim([edges(1) edges(end)]); ylim([0 max([h12.Values, h21.Values])+0.5]); box on; grid on
    set(gca, 'YTick',unique(ceil(get(gca, 'YTick'))), 'XTick',unique(ceil(get(gca, 'XTick'))), 'FontSize',axesLabelSize);
    xlabel('Number of Times Exchanged', 'FontSize',xyLabelSize)
    ylabel('Number of Extracellular Metabolites', 'FontSize',xyLabelSize)
    title([model_title ': Exchanged Metabolites'], 'FontSize',titleSize)
    l_h = legend(headers, 'FontSize',axesLabelSize, 'Location','NorthEast'); set(l_h, 'EdgeColor','w');
    saveas(gcf,[plotFolder2 'exchMets_notInMedium_hist'],'png');
    saveas(gcf,[plotFolder2 'exchMets_notInMedium_hist'],'fig');
end

end

