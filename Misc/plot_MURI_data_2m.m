function plot_MURI_data_2m(model2,homePath,saveFolder,model_title)
%PLOT_MURI_DATA Create all the plots for the MURI data
%
% Meghan Thommes 03/20/2017

%% Biomass Flux

axesLabelSize = 10;
xyLabelSize = 18;
titleSize = 22;
lineWidth = 3;
markerSize = 25;

figure(1);
plot(model2{1}.sparse_con,model2{1}.biomass + model2{2}.biomass,'ro-', 'LineWidth',lineWidth); hold on
plot(model2{1}.sparse_con,model2{1}.biomass,'gd-', 'LineWidth',lineWidth);
plot(model2{2}.sparse_con,model2{2}.biomass,'bs-', 'LineWidth',lineWidth); hold off
set(gca, 'FontSize',axesLabelSize); grid on
xlabel('Constraint Bound for the Number of Intracellular/Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Biomass Flux', 'FontSize',xyLabelSize)
title(model_title, 'FontSize',titleSize)
l_h = legend('Total Biomass of 2 Models','Biomass of Model #1','Biomass of Model #2', 'Location','East'); set(l_h, 'EdgeColor','w');
saveas(gcf,[homePath saveFolder 'biomass_flux'],'png'); saveas(gcf,[homePath saveFolder 'biomass_flux'],'fig');

%% Binary Variables

sparseCon_idx_2m = find(model2{1}.biomass~=0);

figure(3);
imagesc(model2{1}.int(:,sparseCon_idx_2m))
colormap(parula(2)); c = colorbar;
c.Label.String = 'Binary Variable'; c.Label.FontSize = xyLabelSize;
c.Ticks = [0.25 0.75]; c.TickLabels = [0,1]; c.FontSize = axesLabelSize;
set(gca, 'YTick',1:numel(model2{1}.intl_idx), 'YTickLabel',strrep(model2{1}.rxns(model2{1}.intl_idx),'_',' '), ...
    'XTick',1:numel(sparseCon_idx_2m), 'XTickLabel',model2{1}.sparse_con(sparseCon_idx_2m), ...
    'FontSize',10)
xlabel('Constraint Bound for the Number of Intracellular/Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Intracellular/Internal Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 1'], 'FontSize',titleSize)
saveas(gcf,[homePath saveFolder 'int_m1'],'png'); saveas(gcf,[homePath saveFolder 'int_m1'],'fig');

figure(4);
imagesc(model2{2}.int(:,sparseCon_idx_2m))
colormap(parula(2)); c = colorbar;
c.Label.String = 'Binary Variable'; c.Label.FontSize = xyLabelSize;
c.Ticks = [0.25 0.75]; c.TickLabels = [0,1]; c.FontSize = axesLabelSize;
set(gca, 'YTick',1:numel(model2{2}.intl_idx), 'YTickLabel',strrep(model2{2}.rxns(model2{2}.intl_idx),'_',' '), ...
    'XTick',1:numel(sparseCon_idx_2m), 'XTickLabel',model2{2}.sparse_con(sparseCon_idx_2m), ...
    'FontSize',10)
xlabel('Constraint Bound for the Number of Intracellular/Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Intracellular/Internal Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 2'], 'FontSize',titleSize)
saveas(gcf,[homePath saveFolder 'int_m2'],'png'); saveas(gcf,[homePath saveFolder 'int_m2'],'fig');

%% Reaction Fluxes

flux_lim = [min([min(min(model2{1}.flux(model2{1}.intl_idx,sparseCon_idx_2m))), ...
    min(min(model2{2}.flux(model2{2}.intl_idx,sparseCon_idx_2m)))]), ...
    max([max(max(model2{1}.flux(model2{1}.intl_idx,sparseCon_idx_2m))), ...
    max(max(model2{2}.flux(model2{2}.intl_idx,sparseCon_idx_2m)))])];

X = model2{1}.flux(model2{1}.intl_idx,sparseCon_idx_2m);
X(model2{1}.int(:,sparseCon_idx_2m) == 0) = NaN; % set t_i=0 to NaN
figure(6);
h = imagesc(X); set(h,'AlphaData',~isnan(X)); % make NaNs white
caxis(flux_lim); c_mets = colorbar;
c_mets.Label.String = 'Reaction Flux (mmol/gCDW/hr)'; c_mets.Label.FontSize = xyLabelSize;
set(gca, 'YTick',1:numel(model2{1}.intl_idx), 'YTickLabel',strrep(model2{1}.rxns(model2{1}.intl_idx),'_',' '), ...
    'XTick',1:numel(sparseCon_idx_2m), 'XTickLabel',model2{1}.sparse_con(sparseCon_idx_2m), ...
    'FontSize',10)
xlabel('Constraint Bound for the Number of Intracellular/Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Intracellular/Internal Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 1'], 'FontSize',titleSize)
saveas(gcf,[homePath saveFolder 'flux_m1'],'png'); saveas(gcf,[homePath saveFolder 'flux_m1'],'fig');

X = model2{2}.flux(model2{2}.intl_idx,sparseCon_idx_2m);
X(model2{2}.int(:,sparseCon_idx_2m) == 0) = NaN; % set t_i=0 to NaN
figure(7);
h = imagesc(X); set(h,'AlphaData',~isnan(X)); % make NaNs white
caxis(flux_lim); c_mets = colorbar;
c_mets.Label.String = 'Reaction Flux (mmol/gCDW/hr)'; c_mets.Label.FontSize = xyLabelSize;
set(gca, 'YTick',1:numel(model2{2}.intl_idx), 'YTickLabel',strrep(model2{2}.rxns(model2{2}.intl_idx),'_',' '), ...
    'XTick',1:numel(sparseCon_idx_2m), 'XTickLabel',model2{2}.sparse_con(sparseCon_idx_2m), ...
    'FontSize',10)
xlabel('Constraint Bound for the Number of Intracellular/Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Intracellular/Internal Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 2'], 'FontSize',titleSize)
saveas(gcf,[homePath saveFolder 'flux_m2'],'png'); saveas(gcf,[homePath saveFolder 'flux_m2'],'fig');

%% Plot Metabolite Exchange

arrow_scale = 4;
% Find Metabolites with Non-Zero Flux
[plot_idx_medium,~] = find(model2{1}.fluxType ~= 0); plot_idx_medium = (unique(plot_idx_medium));
[~,~,med_idx] = intersect(model2{1}.medium,model2{1}.trspt_mets);
% Find Metabolites with Non-Zero Flux that aren't in the Medium
plot_idx_noMedium = setdiff(plot_idx_medium,med_idx);

c_mets = cbrewer('div','Spectral',numel(model2{1}.trspt_mets));
% Flux Exchange
fig = figure(8);
[~,~,mets_handle1] = plotCocultureExchange(model2{1}.trspt_flux(:,sparseCon_idx_2m)',model2{2}.trspt_flux(:,sparseCon_idx_2m)',1,arrow_scale,plot_idx_medium,fig,c_mets(plot_idx_medium,:));
% xlim(flux_lim); ylim(flux_lim);
set(gca, 'FontSize',axesLabelSize); grid on
xlabel('Model 1 Extracellular Metabolite Flux', 'FontSize',xyLabelSize)
ylabel('Model 2 Extracellular Metabolite Flux', 'FontSize',xyLabelSize)
legend(mets_handle1,model2{1}.trspt_mets(plot_idx_medium), 'Interpreter','none', 'Location','NorthEastOutside', 'Box','Off')
title(['2 ' model_title ' Models: Metabolite Exchange'], 'FontSize',titleSize)
saveas(gcf,[homePath saveFolder 'exchange_flux'],'png'); saveas(gcf,[homePath saveFolder 'exchange_flux'],'fig');

% Flux Exchange - No Metabolites from the Medium
fig = figure(9);
[~,~,mets_handle2] = plotCocultureExchange(model2{1}.trspt_flux(:,sparseCon_idx_2m)',model2{2}.trspt_flux(:,sparseCon_idx_2m)',1,arrow_scale,plot_idx_noMedium,fig,c_mets(plot_idx_noMedium,:));
set(gca, 'FontSize',axesLabelSize); grid on
xlabel('Model 1 Extracellular Metabolite Flux', 'FontSize',xyLabelSize)
ylabel('Model 2 Extracellular Metabolite Flux', 'FontSize',xyLabelSize)
legend(mets_handle2,model2{1}.trspt_mets(plot_idx_noMedium), 'Interpreter','none', 'Location','NorthEastOutside', 'Box','Off')
title(['2 ' model_title ' Models: Metabolite Exchange'], 'FontSize',titleSize)
saveas(gcf,[homePath saveFolder 'exchange_flux_withoutMediumMets'],'png'); saveas(gcf,[homePath saveFolder 'exchange_flux_withoutMediumMets'],'fig');

% LogScale
X = sign(model2{1}.trspt_flux(:,sparseCon_idx_2m)').*log10(abs(model2{1}.trspt_flux(:,sparseCon_idx_2m)')); X(isnan(X)) = 0;
Y = sign(model2{2}.trspt_flux(:,sparseCon_idx_2m)').*log10(abs(model2{2}.trspt_flux(:,sparseCon_idx_2m)')); Y(isnan(Y)) = 0;
fig = figure(10);
[~,~,mets_handle2] = plotCocultureExchange(X,Y,1,1.5,plot_idx_noMedium,fig,c_mets(plot_idx_noMedium,:));
set(gca, 'FontSize',axesLabelSize); grid on
xlabel('Model 1 Extracellular Metabolite Flux (log10)', 'FontSize',xyLabelSize)
ylabel('Model 2 Extracellular Metabolite Flux (log10)', 'FontSize',xyLabelSize)
legend(mets_handle2,model2{1}.trspt_mets(plot_idx_noMedium), 'Interpreter','none', 'Location','NorthEastOutside', 'Box','Off')
title(['2 ' model_title ' Models: Metabolite Exchange'], 'FontSize',titleSize)
saveas(gcf,[homePath saveFolder 'exchange_flux_withoutMediumMets_log10'],'png'); saveas(gcf,[homePath saveFolder 'exchange_flux_withoutMediumMets_log10'],'fig');

%% Plot Extracellular Flux

medium_lim = [floor(min(min([model2{1}.trspt_flux(plot_idx_medium,sparseCon_idx_2m), model2{2}.trspt_flux(plot_idx_medium,sparseCon_idx_2m)]))); ...
    ceil(max(max([model2{1}.trspt_flux(plot_idx_medium,sparseCon_idx_2m), model2{2}.trspt_flux(plot_idx_medium,sparseCon_idx_2m)])));];
noMedium_lim = [floor(min(min([model2{1}.trspt_flux(plot_idx_noMedium,sparseCon_idx_2m), model2{2}.trspt_flux(plot_idx_noMedium,sparseCon_idx_2m)]))); ...
    ceil(max(max([model2{1}.trspt_flux(plot_idx_noMedium,sparseCon_idx_2m), model2{2}.trspt_flux(plot_idx_noMedium,sparseCon_idx_2m)])));];

figure(13);
set(gcf,'DefaultAxesColorOrder',c_mets(plot_idx_medium,:))
plot(repmat(model2{1}.sparse_con(sparseCon_idx_2m),numel(plot_idx_medium),1)',model2{1}.trspt_flux(plot_idx_medium,sparseCon_idx_2m)', 'o-', 'LineWidth',lineWidth)
ylim(medium_lim); grid on; box on
xlabel('Constraint Bound for the Number of Intracellular/Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Extracellular/External Metabolite Flux', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 1 Metabolite Exchange'], 'FontSize',titleSize)
legend(model2{1}.trspt_mets(plot_idx_medium), 'Location','NorthEastOutside', 'Box','Off')
saveas(gcf,[homePath saveFolder 'extracellular_flux_2models_m1'],'png'); saveas(gcf,[homePath saveFolder 'extracellular_flux_2models_m1'],'fig');

figure(14);
set(gcf,'DefaultAxesColorOrder',c_mets(plot_idx_noMedium,:))
plot(repmat(model2{1}.sparse_con(sparseCon_idx_2m),numel(plot_idx_noMedium),1)',model2{1}.trspt_flux(plot_idx_noMedium,sparseCon_idx_2m)', 'o-', 'LineWidth',lineWidth)
ylim(noMedium_lim); grid on; box on
xlabel('Constraint Bound for the Number of Intracellular/Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Extracellular/External Metabolite Flux', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 1 Metabolite Exchange'], 'FontSize',titleSize)
legend(model2{1}.trspt_mets(plot_idx_noMedium), 'Location','NorthEastOutside', 'Box','Off')
saveas(gcf,[homePath saveFolder 'extracellular_flux_2models_m1_withoutMediumMets'],'png'); saveas(gcf,[homePath saveFolder 'extracellular_flux_2models_m1_withoutMediumMets'],'fig');

figure(15);
set(gcf,'DefaultAxesColorOrder',c_mets(plot_idx_medium,:))
plot(repmat(model2{2}.sparse_con(sparseCon_idx_2m),numel(plot_idx_medium),1)',model2{2}.trspt_flux(plot_idx_medium,sparseCon_idx_2m)', 'o-', 'LineWidth',lineWidth)
ylim(medium_lim); grid on; box on
xlabel('Constraint Bound for the Number of Intracellular/Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Extracellular/External Metabolite Flux', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 2 Metabolite Exchange'], 'FontSize',titleSize)
legend(model2{2}.trspt_mets(plot_idx_medium), 'Location','NorthEastOutside', 'Box','Off')
saveas(gcf,[homePath saveFolder 'extracellular_flux_2models_m2'],'png'); saveas(gcf,[homePath saveFolder 'extracellular_flux_2models_m2'],'fig');

figure(16);
set(gcf,'DefaultAxesColorOrder',c_mets(plot_idx_noMedium,:))
plot(repmat(model2{2}.sparse_con(sparseCon_idx_2m),numel(plot_idx_noMedium),1)',model2{2}.trspt_flux(plot_idx_noMedium,sparseCon_idx_2m)', 'o-', 'LineWidth',lineWidth)
ylim(noMedium_lim); grid on; box on
xlabel('Constraint Bound for the Number of Intracellular/Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Extracellular/External Metabolite Flux', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 2 Metabolite Exchange'], 'FontSize',titleSize)
legend(model2{2}.trspt_mets(plot_idx_noMedium), 'Location','NorthEastOutside', 'Box','Off')
saveas(gcf,[homePath saveFolder 'extracellular_flux_2models_m2_withoutMediumMets'],'png'); saveas(gcf,[homePath saveFolder 'extracellular_flux_2models_m2_withoutMediumMets'],'fig');

%% Plot Exchange Type

cmap = [[0.90,0.90,0.90]; ... % 0:light gray (not used by either model)
    [0.00,0.30,0.10]; ... % 1:dark green (produced by both models)
    [0.25,0.55,0.10]; ... % 2:model2{1}.medium green (produced by model 1)
    [0.50,0.80,0.10]; ... % 3:light green (produced by model 2)
    [0.00,0.10,0.30]; ... % 4:dark blue (consumed by both models)
    [0.00,0.35,0.55]; ... % 5:model2{1}.medium blue (consumed by model 1)
    [0.00,0.60,0.80]; ... % 6:light blue (consumed by model 2)
    [1.00,0.80,0.10]; ... % 7:golden yellow (2 -> 1)
    [0.90,0.60,0.10]]; ... % 8:brown yellow (1 -> 2)

% Heat Map of Flux Type: * indicates potential cross-feeding
total_trsptFlux = model2{1}.trspt_flux + model2{2}.trspt_flux;
[nozeros,c] = find((model2{1}.fluxType(:,sparseCon_idx_2m)==7 | model2{1}.fluxType(:,sparseCon_idx_2m)==8) & total_trsptFlux(:,sparseCon_idx_2m)>=0);
R = setdiff(unique(nozeros),med_idx);
C = c(ismember(nozeros,R));
R = nozeros(ismember(nozeros,R));

figure(17);
imagesc(model2{1}.fluxType(:,sparseCon_idx_2m)); colormap(cmap);
c = colorbar; c.Label.String = 'Exchange Type'; c.Label.FontSize = xyLabelSize;
c.Ticks = 0.45:8/9:8; c.TickLabels = 0:1:8;
hold on; plot(C,R,'k*', 'LineWidth',1.1, 'MarkerSize',8); hold off
set(gca, 'YTick',1:numel(model2{1}.trspt_mets), 'YTickLabel',model2{1}.trspt_mets, ...
    'XTick',1:numel(sparseCon_idx_2m), 'XTickLabel',model2{1}.sparse_con(sparseCon_idx_2m))
xlabel('Constraint Bound for the Number of Intracellular/Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Extracellular/External Metabolites', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Metabolite Exchange'], 'FontSize',titleSize)
saveas(gcf,[homePath saveFolder 'metabolite_type'],'png'); saveas(gcf,[homePath saveFolder 'metabolite_type'],'fig');

% Heat Map of Exchange Flux: * indicates potential cross-feeding
figure(18);
imagesc(total_trsptFlux(:,sparseCon_idx_2m));
c = colorbar; c.Label.String = 'Total Metabolite Flux (mmol/gCDW/hr)'; c.Label.FontSize = xyLabelSize;
hold on; plot(C,R,'k*', 'LineWidth',1.1, 'MarkerSize',8); hold off
set(gca, 'YTick',1:numel(model2{1}.trspt_mets), 'YTickLabel',model2{1}.trspt_mets, ...
    'XTick',1:numel(sparseCon_idx_2m), 'XTickLabel',model2{1}.sparse_con(sparseCon_idx_2m))
xlabel('Constraint Bound for the Number of Intracellular/Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Extracellular/External Metabolites', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Extracellular Metabolite Flux'], 'FontSize',titleSize)
saveas(gcf,[homePath saveFolder 'metabolite_totalFlux'],'png'); saveas(gcf,[homePath saveFolder 'metabolite_totalFlux'],'fig');

%% Cluster - 2 Models - Rxns

sparseCon_idx_2m_noUnconstrained = sparseCon_idx_2m(1:end-1);

% Int - Biclustered by Binary Variable, All Reactions
int_data = [model2{1}.int(:,sparseCon_idx_2m_noUnconstrained); model2{2}.int(:,sparseCon_idx_2m_noUnconstrained)];
xnames = cellstr(int2str(model2{1}.sparse_con(sparseCon_idx_2m_noUnconstrained)'));
ynames1 = arrayfun(@(x) [model2{1}.rxns{model2{1}.intl_idx(x)} '-1'],1:numel(model2{1}.intl_idx),'Uni',false);
ynames2 = arrayfun(@(x) [model2{2}.rxns{model2{2}.intl_idx(x)} '-2'],1:numel(model2{2}.intl_idx),'Uni',false);
ynames = strrep([ynames1, ynames2],'_',' ');
fig_handle(1) = figure(25); fig_handle(2) = figure(26);
[P1,P2,axes_handle] = cluster_flux(int_data,xnames,ynames,fig_handle);
% Dendrogram
set(axes_handle(1), 'FontSize',6); set(get(axes_handle(1), 'YLabel'), 'FontSize',xyLabelSize)
xlabel(axes_handle(1),'Reactions', 'FontSize',xyLabelSize)
xlabel(axes_handle(2),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
title(axes_handle(1),['2 ' model_title ' Models: Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
% Heat Map
set(axes_handle(3), 'FontSize',6)
xlabel(axes_handle(3),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel(axes_handle(3),'Reactions', 'FontSize',xyLabelSize)
title(axes_handle(3),['2 ' model_title ' Models: Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
colormap(parula(2)); c = colorbar;
c.Label.String = 'Binary Variable'; c.Label.FontSize = xyLabelSize;
c.Ticks = [0.25 0.75]; c.TickLabels = [0,1]; c.FontSize = axesLabelSize;
saveas(fig_handle(1),[homePath saveFolder 'dendrogram_int_2models_rxns'],'png'); saveas(gcf,[homePath saveFolder 'dendrogram_int_2models_rxns'],'fig');
saveas(fig_handle(2),[homePath saveFolder 'bicluster_int_2models_rxns'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_int_2models_rxns'],'fig');

% Flux - Biclustered by Binary Variable, All Reactions
flux_data = [model2{1}.flux(model2{1}.intl_idx,sparseCon_idx_2m_noUnconstrained); model2{2}.flux(model2{2}.intl_idx,sparseCon_idx_2m_noUnconstrained)];
flux_data(int_data == 0) = NaN; % set t_i=0 to NaN
figure(27)
h = imagesc(flux_data(P1,P2));
set(h,'AlphaData',~isnan(flux_data(P1,P2))); % make NaNs white
set(gca, 'XTick',1:numel(P2), 'XTickLabel',xnames(P2), ...
    'YTick',1:numel(P1), 'YTickLabel',ynames(P1), 'Fontsize',6);
xlabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
caxis(flux_lim); colormap(parula); c = colorbar; c.FontSize = axesLabelSize;
c.Label.String = 'Reaction Flux (mmol/gCDW/hr)'; c.Label.FontSize = xyLabelSize;
saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_rxns'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_rxns'],'fig');

% Int - Biclustered by Binary Variable, No All-Zero Reactions
new_int_data = int_data;
tol = 1E-6; new_int_data(new_int_data>-tol & new_int_data<tol) = 0;
[nozeros,~] = find(new_int_data~=0); nozeros = unique(nozeros);
new_int_data = new_int_data(nozeros,:); new_ynames = ynames(nozeros);
fig_handle(1) = figure(28); fig_handle(2) = figure(29);
[P1,P2,axes_handle] = cluster_flux(new_int_data,xnames,new_ynames,fig_handle);
% Dendrogram
set(axes_handle(1), 'FontSize',10); set(get(axes_handle(1), 'YLabel'), 'FontSize',xyLabelSize)
xlabel(axes_handle(1),'Reactions', 'FontSize',xyLabelSize)
xlabel(axes_handle(2),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
title(axes_handle(1),['2 ' model_title ' Models: Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
% Heat Map
set(axes_handle(3), 'FontSize',10)
xlabel(axes_handle(3),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel(axes_handle(3),'Reactions', 'FontSize',xyLabelSize)
title(axes_handle(3),['2 ' model_title ' Models: Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
colormap(parula(2)); c = colorbar;
c.Label.String = 'Binary Variable'; c.Label.FontSize = xyLabelSize;
c.Ticks = [0.25 0.75]; c.TickLabels = [0,1]; c.FontSize = axesLabelSize;
saveas(fig_handle(1),[homePath saveFolder 'dendrogram_int_2models_rxns_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'dendrogram_int_2models_rxns_nozeros'],'fig');
saveas(fig_handle(2),[homePath saveFolder 'bicluster_int_2models_rxns_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_int_2models_rxns_nozeros'],'fig');

% Flux - Biclustered by Binary Variable, All Intracellular Reactions
figure(30)
new_flux_data = flux_data(nozeros,:);
h = imagesc(new_flux_data(P1,P2));
set(h,'AlphaData',~isnan(new_flux_data(P1,P2))); % make NaNs white
set(gca, 'XTick',1:numel(P2), 'XTickLabel',xnames(P2), ...
    'YTick',1:numel(P1), 'YTickLabel',ynames(P1), 'Fontsize',10);
xlabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
caxis(flux_lim); colormap(parula); c = colorbar; c.FontSize = axesLabelSize;
c.Label.String = 'Reaction Flux (mmol/gCDW/hr)'; c.Label.FontSize = xyLabelSize;
saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_rxns_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_rxns_nozeros'],'fig');

%% Cluster - 2 Models, Model 1 - Rxns

% Int - Biclustered by Binary Variable, All Reactions
int_data = model2{1}.int(:,sparseCon_idx_2m_noUnconstrained);
xnames = cellstr(int2str(model2{1}.sparse_con(sparseCon_idx_2m_noUnconstrained)'));
ynames = strrep(model2{1}.rxns(model2{2}.intl_idx),'_',' ');
fig_handle(1) = figure(31); fig_handle(2) = figure(32);
[P1,P2,axes_handle] = cluster_flux(int_data,xnames,ynames,fig_handle);
% Dendrogram
set(axes_handle(1), 'FontSize',6); set(get(axes_handle(1), 'YLabel'), 'FontSize',xyLabelSize)
xlabel(axes_handle(1),'Reactions', 'FontSize',xyLabelSize)
xlabel(axes_handle(2),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
title(axes_handle(1),['2 ' model_title ' Models: Model 1 Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
% Heat Map
set(axes_handle(3), 'FontSize',6)
xlabel(axes_handle(3),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel(axes_handle(3),'Reactions', 'FontSize',xyLabelSize)
title(axes_handle(3),['2 ' model_title ' Models: Model 1 Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
colormap(parula(2)); c = colorbar;
c.Label.String = 'Binary Variable'; c.Label.FontSize = xyLabelSize;
c.Ticks = [0.25 0.75]; c.TickLabels = [0,1]; c.FontSize = axesLabelSize;
saveas(fig_handle(1),[homePath saveFolder 'dendrogram_int_2models_m1_rxns'],'png'); saveas(gcf,[homePath saveFolder 'dendrogram_int_2models_m1_rxns'],'fig');
saveas(fig_handle(2),[homePath saveFolder 'bicluster_int_2models_m1_rxns'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_int_2models_m1_rxns'],'fig');

% Flux - Biclustered by Binary Variable, All Reactions
flux_data = [model2{1}.flux(model2{1}.intl_idx,sparseCon_idx_2m_noUnconstrained); model2{2}.flux(model2{2}.intl_idx,sparseCon_idx_2m_noUnconstrained)];
flux_data(int_data == 0) = NaN; % set t_i=0 to NaN
figure(33)
h = imagesc(flux_data(P1,P2));
set(h,'AlphaData',~isnan(flux_data(P1,P2))); % make NaNs white
set(gca, 'XTick',1:numel(P2), 'XTickLabel',xnames(P2), ...
    'YTick',1:numel(P1), 'YTickLabel',ynames(P1), 'Fontsize',6);
xlabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 1 Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
caxis(flux_lim); colormap(parula); c = colorbar; c.FontSize = axesLabelSize;
c.Label.String = 'Reaction Flux (mmol/gCDW/hr)'; c.Label.FontSize = xyLabelSize;
saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_m1_rxns'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_m1_rxns'],'fig');

% Int - Biclustered by Binary Variable, No All-Zero Reactions
new_int_data = int_data;
tol = 1E-6; new_int_data(new_int_data>-tol & new_int_data<tol) = 0;
[nozeros,~] = find(new_int_data~=0); nozeros = unique(nozeros);
new_int_data = new_int_data(nozeros,:); new_ynames = ynames(nozeros);
fig_handle(1) = figure(34); fig_handle(2) = figure(35);
[P1,P2,axes_handle] = cluster_flux(new_int_data,xnames,new_ynames,fig_handle);
% Dendrogram
set(axes_handle(1), 'FontSize',10); set(get(axes_handle(1), 'YLabel'), 'FontSize',xyLabelSize)
xlabel(axes_handle(1),'Reactions', 'FontSize',xyLabelSize)
xlabel(axes_handle(2),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
title(axes_handle(1),['2 ' model_title ' Models: Model 1 Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
% Heat Map
set(axes_handle(3), 'FontSize',10)
xlabel(axes_handle(3),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel(axes_handle(3),'Reactions', 'FontSize',xyLabelSize)
title(axes_handle(3),['2 ' model_title ' Models: Model 1 Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
colormap(parula(2)); c = colorbar;
c.Label.String = 'Binary Variable'; c.Label.FontSize = xyLabelSize;
c.Ticks = [0.25 0.75]; c.TickLabels = [0,1]; c.FontSize = axesLabelSize;
saveas(fig_handle(1),[homePath saveFolder 'dendrogram_int_2models_m1_rxns_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'dendrogram_int_2models_m1_rxns_nozeros'],'fig');
saveas(fig_handle(2),[homePath saveFolder 'bicluster_int_2models_m1_rxns_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_int_2models_m1_rxns_nozeros'],'fig');

% Flux - Biclustered by Binary Variable, All Intracellular Reactions
figure(36)
new_flux_data = flux_data(nozeros,:);
h = imagesc(new_flux_data(P1,P2));
set(h,'AlphaData',~isnan(new_flux_data(P1,P2))); % make NaNs white
set(gca, 'XTick',1:numel(P2), 'XTickLabel',xnames(P2), ...
    'YTick',1:numel(P1), 'YTickLabel',ynames(P1), 'Fontsize',10);
xlabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 1 Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
caxis(flux_lim); colormap(parula); c = colorbar; c.FontSize = axesLabelSize;
c.Label.String = 'Reaction Flux (mmol/gCDW/hr)'; c.Label.FontSize = xyLabelSize;
saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_m1_rxns_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_m1_rxns_nozeros'],'fig');

%% Cluster - 2 Models, Model 2 - Rxns

% Int - Biclustered by Binary Variable, All Reactions
int_data = model2{2}.int(:,sparseCon_idx_2m_noUnconstrained);
xnames = cellstr(int2str(model2{2}.sparse_con(sparseCon_idx_2m_noUnconstrained)'));
ynames = strrep(model2{2}.rxns(model2{2}.intl_idx),'_',' ');
fig_handle(1) = figure(37); fig_handle(2) = figure(38);
[P1,P2,axes_handle] = cluster_flux(int_data,xnames,ynames,fig_handle);
% Dendrogram
set(axes_handle(1), 'FontSize',6); set(get(axes_handle(1), 'YLabel'), 'FontSize',xyLabelSize)
xlabel(axes_handle(1),'Reactions', 'FontSize',xyLabelSize)
xlabel(axes_handle(2),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
title(axes_handle(1),['2 ' model_title ' Models: Model 2 Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
% Heat Map
set(axes_handle(3), 'FontSize',6)
xlabel(axes_handle(3),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel(axes_handle(3),'Reactions', 'FontSize',xyLabelSize)
title(axes_handle(3),['2 ' model_title ' Models: Model 2 Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
colormap(parula(2)); c = colorbar;
c.Label.String = 'Binary Variable'; c.Label.FontSize = xyLabelSize;
c.Ticks = [0.25 0.75]; c.TickLabels = [0,1]; c.FontSize = axesLabelSize;
saveas(fig_handle(1),[homePath saveFolder 'dendrogram_int_2models_m2_rxns'],'png'); saveas(gcf,[homePath saveFolder 'dendrogram_int_2models_m2_rxns'],'fig');
saveas(fig_handle(2),[homePath saveFolder 'bicluster_int_2models_m2_rxns'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_int_2models_m2_rxns'],'fig');

% Flux - Biclustered by Binary Variable, All Reactions
flux_data = [model2{1}.flux(model2{1}.intl_idx,sparseCon_idx_2m_noUnconstrained); model2{2}.flux(model2{2}.intl_idx,sparseCon_idx_2m_noUnconstrained)];
flux_data(int_data == 0) = NaN; % set t_i=0 to NaN
figure(39)
h = imagesc(flux_data(P1,P2));
set(h,'AlphaData',~isnan(flux_data(P1,P2))); % make NaNs white
set(gca, 'XTick',1:numel(P2), 'XTickLabel',xnames(P2), ...
    'YTick',1:numel(P1), 'YTickLabel',ynames(P1), 'Fontsize',6);
xlabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 2 Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
caxis(flux_lim); colormap(parula); c = colorbar; c.FontSize = axesLabelSize;
c.Label.String = 'Reaction Flux (mmol/gCDW/hr)'; c.Label.FontSize = xyLabelSize;
saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_m2_rxns'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_m2_rxns'],'fig');

% Int - Biclustered by Binary Variable, No All-Zero Reactions
new_int_data = int_data;
tol = 1E-6; new_int_data(new_int_data>-tol & new_int_data<tol) = 0;
[nozeros,~] = find(new_int_data~=0); nozeros = unique(nozeros);
new_int_data = new_int_data(nozeros,:); new_ynames = ynames(nozeros);
fig_handle(1) = figure(40); fig_handle(2) = figure(41);
[P1,P2,axes_handle] = cluster_flux(new_int_data,xnames,new_ynames,fig_handle);
% Dendrogram
set(axes_handle(1), 'FontSize',10); set(get(axes_handle(1), 'YLabel'), 'FontSize',xyLabelSize)
xlabel(axes_handle(1),'Reactions', 'FontSize',xyLabelSize)
xlabel(axes_handle(2),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
title(axes_handle(1),['2 ' model_title ' Models: Model 2 Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
% Heat Map
set(axes_handle(3), 'FontSize',10)
xlabel(axes_handle(3),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel(axes_handle(3),'Reactions', 'FontSize',xyLabelSize)
title(axes_handle(3),['2 ' model_title ' Models: Model 2 Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
colormap(parula(2)); c = colorbar;
c.Label.String = 'Binary Variable'; c.Label.FontSize = xyLabelSize;
c.Ticks = [0.25 0.75]; c.TickLabels = [0,1]; c.FontSize = axesLabelSize;
saveas(fig_handle(1),[homePath saveFolder 'dendrogram_int_2models_m2_rxns_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'dendrogram_int_2models_m2_rxns_nozeros'],'fig');
saveas(fig_handle(2),[homePath saveFolder 'bicluster_int_2models_m2_rxns_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_int_2models_m2_rxns_nozeros'],'fig');

% Flux - Biclustered by Binary Variable, All Intracellular Reactions
figure(42)
new_flux_data = flux_data(nozeros,:);
h = imagesc(new_flux_data(P1,P2));
set(h,'AlphaData',~isnan(new_flux_data(P1,P2))); % make NaNs white
set(gca, 'XTick',1:numel(P2), 'XTickLabel',xnames(P2), ...
    'YTick',1:numel(P1), 'YTickLabel',ynames(P1), 'Fontsize',10);
xlabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 2 Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
caxis(flux_lim); colormap(parula); c = colorbar; c.FontSize = axesLabelSize;
c.Label.String = 'Reaction Flux (mmol/gCDW/hr)'; c.Label.FontSize = xyLabelSize;
saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_m2_rxns_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_m2_rxns_nozeros']);

%% Cluster - 2 Models - SparseCon

% Int - Biclustered by Binary Variable, All Reactions
int_data = [model2{1}.int(:,sparseCon_idx_2m_noUnconstrained), model2{2}.int(:,sparseCon_idx_2m_noUnconstrained)];
xnames1 = cellstr(int2str(model2{1}.sparse_con(sparseCon_idx_2m_noUnconstrained)'));
xnames2 = cellstr(int2str(model2{2}.sparse_con(sparseCon_idx_2m_noUnconstrained)'));
xnames1 = arrayfun(@(x) [xnames1{x} '-1'],1:numel(xnames1),'Uni',false);
xnames2 = arrayfun(@(x) [xnames2{x} '-2'],1:numel(xnames2),'Uni',false);
xnames = [xnames1, xnames2];
ynames = strrep(model2{1}.rxns(model2{1}.intl_idx),'_',' ');
fig_handle(1) = figure(43); fig_handle(2) = figure(44);
[P1,P2,axes_handle] = cluster_flux(int_data,xnames,ynames,fig_handle);
% Dendrogram
set(axes_handle(1), 'FontSize',10); set(get(axes_handle(1), 'YLabel'), 'FontSize',xyLabelSize)
xlabel(axes_handle(1),'Reactions', 'FontSize',0.75*xyLabelSize)
xlabel(axes_handle(2),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
title(axes_handle(1),['2 ' model_title ' Models: Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
% Heat Map
set(axes_handle(3), 'FontSize',10)
xlabel(axes_handle(3),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel(axes_handle(3),'Reactions', 'FontSize',xyLabelSize)
title(axes_handle(3),['2 ' model_title ' Models: Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
colormap(parula(2)); c = colorbar;
c.Label.String = 'Binary Variable'; c.Label.FontSize = xyLabelSize;
c.Ticks = [0.25 0.75]; c.TickLabels = [0,1]; c.FontSize = axesLabelSize;
saveas(fig_handle(1),[homePath saveFolder 'dendrogram_int_2models_sparseCon'],'png'); saveas(gcf,[homePath saveFolder 'dendrogram_int_2models_sparseCon'],'fig');
saveas(fig_handle(2),[homePath saveFolder 'bicluster_int_2models_sparseCon'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_int_2models_sparseCon'],'fig');

% Flux - Biclustered by Binary Variable, All Reactions
flux_data = [model2{1}.flux(model2{1}.intl_idx,sparseCon_idx_2m_noUnconstrained), model2{2}.flux(model2{2}.intl_idx,sparseCon_idx_2m_noUnconstrained)];
flux_data(int_data == 0) = NaN; % set t_i=0 to NaN
figure(45)
h = imagesc(flux_data(P1,P2));
set(h,'AlphaData',~isnan(flux_data(P1,P2))); % make NaNs white
set(gca, 'XTick',1:numel(P2), 'XTickLabel',xnames(P2), ...
    'YTick',1:numel(P1), 'YTickLabel',ynames(P1), 'Fontsize',10);
xlabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
c = colorbar; c.Label.String = 'Reaction Flux (mmol/gCDW/hr)';
c.Label.FontSize = xyLabelSize; c.FontSize = axesLabelSize;
saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_sparseCon'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_sparseCon'],'fig');

% Int - Biclustered by Binary Variable, No All-Zero Reactions
new_int_data = int_data;
tol = 1E-6; new_int_data(new_int_data>-tol & new_int_data<tol) = 0;
[nozeros,~] = find(new_int_data~=0); nozeros = unique(nozeros);
new_int_data = new_int_data(nozeros,:); new_ynames = ynames(nozeros);
fig_handle(1) = figure(46); fig_handle(2) = figure(47);
[P1,P2,axes_handle] = cluster_flux(new_int_data,xnames,new_ynames,fig_handle);
% Dendrogram
set(axes_handle(1), 'FontSize',10); set(get(axes_handle(1), 'YLabel'), 'FontSize',xyLabelSize)
xlabel(axes_handle(1),'Reactions', 'FontSize',xyLabelSize)
xlabel(axes_handle(2),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
title(axes_handle(1),['2 ' model_title ' Models: Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
% Heat Map
set(axes_handle(3), 'FontSize',10)
xlabel(axes_handle(3),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel(axes_handle(3),'Reactions', 'FontSize',xyLabelSize)
title(axes_handle(3),['2 ' model_title ' Models: Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
colormap(parula(2)); c = colorbar;
c.Label.String = 'Binary Variable'; c.Label.FontSize = xyLabelSize;
c.Ticks = [0.25 0.75]; c.TickLabels = [0,1]; c.FontSize = axesLabelSize;
saveas(fig_handle(1),[homePath saveFolder 'dendrogram_int_2models_sparseCon_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'dendrogram_int_2models_sparseCon_nozeros'],'fig');
saveas(fig_handle(2),[homePath saveFolder 'bicluster_int_2models_sparseCon_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_int_2models_sparseCon_nozeros'],'fig');

% Flux - Biclustered by Binary Variable, No All-Zero Reactions
figure(48);
new_flux_data = flux_data(nozeros,:);
h = imagesc(new_flux_data(P1,P2));
set(h,'AlphaData',~isnan(new_flux_data(P1,P2))); % make NaNs white
set(gca, 'XTick',1:numel(P2), 'XTickLabel',xnames(P2), ...
    'YTick',1:numel(P1), 'YTickLabel',ynames(P1), 'Fontsize',10);
xlabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
caxis(flux_lim); colormap(parula); c = colorbar; c.FontSize = axesLabelSize;
c.Label.String = 'Reaction Flux (mmol/gCDW/hr)'; c.Label.FontSize = xyLabelSize;
saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_sparseCon_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_sparseCon_nozeros'],'fig');

%% Cluster - 2 Models, Model 1 - SparseCon

% Int - Biclustered by Binary Variable, All Reactions
int_data = model2{1}.int(:,sparseCon_idx_2m_noUnconstrained);
xnames = cellstr(int2str(model2{1}.sparse_con(sparseCon_idx_2m_noUnconstrained)'));
ynames = strrep(model2{1}.rxns(model2{1}.intl_idx),'_',' ');
fig_handle(1) = figure(49); fig_handle(2) = figure(50);
[P1,P2,axes_handle] = cluster_flux(int_data,xnames,ynames,fig_handle);
% Dendrogram
set(axes_handle(1), 'FontSize',10); set(get(axes_handle(1), 'YLabel'), 'FontSize',xyLabelSize)
xlabel(axes_handle(1),'Reactions', 'FontSize',0.75*xyLabelSize)
xlabel(axes_handle(2),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
title(axes_handle(1),['2 ' model_title ' Models: Model 1 Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
% Heat Map
set(axes_handle(3), 'FontSize',10)
xlabel(axes_handle(3),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel(axes_handle(3),'Reactions', 'FontSize',xyLabelSize)
title(axes_handle(3),['2 ' model_title ' Models: Model 1 Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
colormap(parula(2)); c = colorbar;
c.Label.String = 'Binary Variable'; c.Label.FontSize = xyLabelSize;
c.Ticks = [0.25 0.75]; c.TickLabels = [0,1]; c.FontSize = axesLabelSize;
saveas(fig_handle(1),[homePath saveFolder 'dendrogram_int_2models_m1_sparseCon'],'png'); saveas(gcf,[homePath saveFolder 'dendrogram_int_2models_m1_sparseCon'],'fig');
saveas(fig_handle(2),[homePath saveFolder 'bicluster_int_2models_m1_sparseCon'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_int_2models_m1_sparseCon'],'fig');

% Flux - Biclustered by Binary Variable, All Reactions
flux_data = model2{1}.flux(model2{1}.intl_idx,sparseCon_idx_2m_noUnconstrained);
flux_data(int_data == 0) = NaN; % set t_i=0 to NaN
figure(51)
h = imagesc(flux_data(P1,P2));
set(h,'AlphaData',~isnan(flux_data(P1,P2))); % make NaNs white
set(gca, 'XTick',1:numel(P2), 'XTickLabel',xnames(P2), ...
    'YTick',1:numel(P1), 'YTickLabel',ynames(P1), 'Fontsize',10);
xlabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 1 Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
c = colorbar; c.Label.String = 'Reaction Flux (mmol/gCDW/hr)';
c.Label.FontSize = xyLabelSize; c.FontSize = axesLabelSize;
saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_m1_sparseCon'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_m1_sparseCon'],'fig');

% Int - Biclustered by Binary Variable, No All-Zero Reactions
new_int_data = int_data;
tol = 1E-6; new_int_data(new_int_data>-tol & new_int_data<tol) = 0;
[nozeros,~] = find(new_int_data~=0); nozeros = unique(nozeros);
new_int_data = new_int_data(nozeros,:); new_ynames = ynames(nozeros);
fig_handle(1) = figure(52); fig_handle(2) = figure(53);
[P1,P2,axes_handle] = cluster_flux(new_int_data,xnames,new_ynames,fig_handle);
% Dendrogram
set(axes_handle(1), 'FontSize',10); set(get(axes_handle(1), 'YLabel'), 'FontSize',xyLabelSize)
xlabel(axes_handle(1),'Reactions', 'FontSize',xyLabelSize)
xlabel(axes_handle(2),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
title(axes_handle(1),['2 ' model_title ' Models: Model 1 Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
% Heat Map
set(axes_handle(3), 'FontSize',10)
xlabel(axes_handle(3),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel(axes_handle(3),'Reactions', 'FontSize',xyLabelSize)
title(axes_handle(3),['2 ' model_title ' Models: Model 1 Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
colormap(parula(2)); c = colorbar;
c.Label.String = 'Binary Variable'; c.Label.FontSize = xyLabelSize;
c.Ticks = [0.25 0.75]; c.TickLabels = [0,1]; c.FontSize = axesLabelSize;
saveas(fig_handle(1),[homePath saveFolder 'dendrogram_int_2models_m1_sparseCon_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'dendrogram_int_2models_m1_sparseCon_nozeros'],'fig');
saveas(fig_handle(2),[homePath saveFolder 'bicluster_int_2models_m1_sparseCon_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_int_2models_m1_sparseCon_nozeros'],'fig');

% Flux - Biclustered by Binary Variable, No All-Zero Reactions
figure(54);
new_flux_data = flux_data(nozeros,:);
h = imagesc(new_flux_data(P1,P2));
set(h,'AlphaData',~isnan(new_flux_data(P1,P2))); % make NaNs white
set(gca, 'XTick',1:numel(P2), 'XTickLabel',xnames(P2), ...
    'YTick',1:numel(P1), 'YTickLabel',ynames(P1), 'Fontsize',10);
xlabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 1 Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
caxis(flux_lim); colormap(parula); c = colorbar; c.FontSize = axesLabelSize;
c.Label.String = 'Reaction Flux (mmol/gCDW/hr)'; c.Label.FontSize = xyLabelSize;
saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_m1_sparseCon_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_m1_sparseCon_nozeros'],'fig');

%% Cluster - 2 Models, Model 2 - SparseCon

% Int - Biclustered by Binary Variable, All Reactions
int_data = model2{2}.int(:,sparseCon_idx_2m_noUnconstrained);
xnames = cellstr(int2str(model2{2}.sparse_con(sparseCon_idx_2m_noUnconstrained)'));
ynames = strrep(model2{2}.rxns(model2{2}.intl_idx),'_',' ');
fig_handle(1) = figure(55); fig_handle(2) = figure(56);
[P1,P2,axes_handle] = cluster_flux(int_data,xnames,ynames,fig_handle);
% Dendrogram
set(axes_handle(1), 'FontSize',10); set(get(axes_handle(1), 'YLabel'), 'FontSize',xyLabelSize)
xlabel(axes_handle(1),'Reactions', 'FontSize',0.75*xyLabelSize)
xlabel(axes_handle(2),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
title(axes_handle(1),['2 ' model_title ' Models: Model 2 Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
% Heat Map
set(axes_handle(3), 'FontSize',10)
xlabel(axes_handle(3),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel(axes_handle(3),'Reactions', 'FontSize',xyLabelSize)
title(axes_handle(3),['2 ' model_title ' Models: Model 2 Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
colormap(parula(2)); c = colorbar;
c.Label.String = 'Binary Variable'; c.Label.FontSize = xyLabelSize;
c.Ticks = [0.25 0.75]; c.TickLabels = [0,1]; c.FontSize = axesLabelSize;
saveas(fig_handle(1),[homePath saveFolder 'dendrogram_int_2models_m2_sparseCon'],'png'); saveas(gcf,[homePath saveFolder 'dendrogram_int_2models_m2_sparseCon'],'fig');
saveas(fig_handle(2),[homePath saveFolder 'bicluster_int_2models_m2_sparseCon'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_int_2models_m2_sparseCon'],'fig');

% Flux - Biclustered by Binary Variable, All Reactions
flux_data = model2{2}.flux(model2{2}.intl_idx,sparseCon_idx_2m_noUnconstrained);
flux_data(int_data == 0) = NaN; % set t_i=0 to NaN
figure(57)
h = imagesc(flux_data(P1,P2));
set(h,'AlphaData',~isnan(flux_data(P1,P2))); % make NaNs white
set(gca, 'XTick',1:numel(P2), 'XTickLabel',xnames(P2), ...
    'YTick',1:numel(P1), 'YTickLabel',ynames(P1), 'Fontsize',10);
xlabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 2 Biclustered by Binary Variable, All Intracellular Reactions'], 'FontSize',titleSize)
c = colorbar; c.Label.String = 'Reaction Flux (mmol/gCDW/hr)';
c.Label.FontSize = xyLabelSize; c.FontSize = axesLabelSize;
saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_m2_sparseCon'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_m2_sparseCon'],'fig');

% Int - Biclustered by Binary Variable, No All-Zero Reactions
new_int_data = int_data;
tol = 1E-6; new_int_data(new_int_data>-tol & new_int_data<tol) = 0;
[nozeros,~] = find(new_int_data~=0); nozeros = unique(nozeros);
new_int_data = new_int_data(nozeros,:); new_ynames = ynames(nozeros);
fig_handle(1) = figure(58); fig_handle(2) = figure(59);
[P1,P2,axes_handle] = cluster_flux(new_int_data,xnames,new_ynames,fig_handle);
% Dendrogram
set(axes_handle(1), 'FontSize',10); set(get(axes_handle(1), 'YLabel'), 'FontSize',xyLabelSize)
xlabel(axes_handle(1),'Reactions', 'FontSize',xyLabelSize)
xlabel(axes_handle(2),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
title(axes_handle(1),['2 ' model_title ' Models: Model 2 Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
% Heat Map
set(axes_handle(3), 'FontSize',10)
xlabel(axes_handle(3),'Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel(axes_handle(3),'Reactions', 'FontSize',xyLabelSize)
title(axes_handle(3),['2 ' model_title ' Models: Model 2 Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
colormap(parula(2)); c = colorbar;
c.Label.String = 'Binary Variable'; c.Label.FontSize = xyLabelSize;
c.Ticks = [0.25 0.75]; c.TickLabels = [0,1]; c.FontSize = axesLabelSize;
saveas(fig_handle(1),[homePath saveFolder 'dendrogram_int_2models_m2_sparseCon_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'dendrogram_int_2models_m2_sparseCon_nozeros'],'fig');
saveas(fig_handle(2),[homePath saveFolder 'bicluster_int_2models_m2_sparseCon_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_int_2models_m2_sparseCon_nozeros'],'fig')

% Flux - Biclustered by Binary Variable, No All-Zero Reactions
figure(60);
new_flux_data = flux_data(nozeros,:);
h = imagesc(new_flux_data(P1,P2));
set(h,'AlphaData',~isnan(new_flux_data(P1,P2))); % make NaNs white
set(gca, 'XTick',1:numel(P2), 'XTickLabel',xnames(P2), ...
    'YTick',1:numel(P1), 'YTickLabel',ynames(P1), 'Fontsize',10);
xlabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 2 Biclustered by Binary Variable, No All-Zero Reactions'], 'FontSize',titleSize)
caxis(flux_lim); colormap(parula); c = colorbar; c.FontSize = axesLabelSize;
c.Label.String = 'Reaction Flux (mmol/gCDW/hr)'; c.Label.FontSize = xyLabelSize;
saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_m2_sparseCon_nozeros'],'png'); saveas(gcf,[homePath saveFolder 'bicluster_flux_2models_m2_sparseCon_nozeros'],'fig');

%% Correlations - 2 Models

% Int
c_int2_sparseCon = corr(model2{1}.int(:,sparseCon_idx_2m_noUnconstrained),model2{2}.int(:,sparseCon_idx_2m_noUnconstrained));
c_int2_rxns = corr(model2{1}.int(:,sparseCon_idx_2m_noUnconstrained)',model2{2}.int(:,sparseCon_idx_2m_noUnconstrained)');

figure(65)
imagesc(c_int2_sparseCon)
set(gca, 'XTick',1:numel(model2{2}.sparse_con(sparseCon_idx_2m_noUnconstrained)), 'XTickLabel',cellstr(int2str(model2{2}.sparse_con(sparseCon_idx_2m_noUnconstrained)')), ...
    'YTick',1:numel(model2{1}.sparse_con(sparseCon_idx_2m_noUnconstrained)), 'YTickLabel',cellstr(int2str(model2{1}.sparse_con(sparseCon_idx_2m_noUnconstrained)')))
xlabel('Constraint Bound for the Number of Internal Reactions (Model 2)', 'FontSize',xyLabelSize)
ylabel('Constraint Bound for the Number of Internal Reactions (Model 1)', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models (by int)'], 'FontSize',titleSize)
c = colorbar; c.Label.String = 'Correlation';
c.Label.FontSize = xyLabelSize; c.FontSize = axesLabelSize;
saveas(gcf,[homePath saveFolder 'corr_int_2models_sparseCon'],'png'); saveas(gcf,[homePath saveFolder 'corr_int_2models_sparseCon'],'fig');

figure(66)
imagesc(c_int2_rxns);
set(gca, 'XTick',1:numel(model2{2}.intl_idx), 'XTickLabel',strrep(model2{2}.rxns(model2{2}.intl_idx),'_',' '), 'XTickLabelRotation',45, ...
    'YTick',1:numel(model2{1}.intl_idx), 'YTickLabel',strrep(model2{1}.rxns(model2{1}.intl_idx),'_',' '), 'FontSize',10)
xlabel('Reactions (Model 2)', 'FontSize',xyLabelSize);
ylabel('Reactions (Model 1)', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models (by int)'], 'FontSize',titleSize)
c = colorbar; c.Label.String = 'Correlation';
c.Label.FontSize = xyLabelSize; c.FontSize = axesLabelSize;
saveas(gcf,[homePath saveFolder 'corr_int_2models_rxns'],'png'); saveas(gcf,[homePath saveFolder 'corr_int_2models_rxns'],'fig');

% Int - Model 1
c_int2m1_sparseCon = corr(model2{1}.int(:,sparseCon_idx_2m_noUnconstrained));
c_int2m1_rxns = corr(model2{1}.int(:,sparseCon_idx_2m_noUnconstrained)');

figure(67)
imagesc(c_int2m1_sparseCon)
set(gca, 'XTick',1:numel(model2{1}.sparse_con(sparseCon_idx_2m_noUnconstrained)), 'XTickLabel',cellstr(int2str(model2{1}.sparse_con(sparseCon_idx_2m_noUnconstrained)')), ...
    'YTick',1:numel(model2{1}.sparse_con(sparseCon_idx_2m_noUnconstrained)), 'YTickLabel',cellstr(int2str(model2{1}.sparse_con(sparseCon_idx_2m_noUnconstrained)')))
xlabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 1 (by int)'], 'FontSize',titleSize)
c = colorbar; c.Label.String = 'Correlation';
c.Label.FontSize = xyLabelSize; c.FontSize = axesLabelSize;
saveas(gcf,[homePath saveFolder 'corr_int_2models_model1_sparseCon'],'png'); saveas(gcf,[homePath saveFolder 'corr_int_2models_model1_sparseCon'],'fig');

figure(68)
imagesc(c_int2m1_rxns);
set(gca, 'XTick',1:numel(model2{1}.intl_idx), 'XTickLabel',strrep(model2{1}.rxns(model2{1}.intl_idx),'_',' '), 'XTickLabelRotation',45, ...
    'YTick',1:numel(model2{1}.intl_idx), 'YTickLabel',strrep(model2{1}.rxns(model2{1}.intl_idx),'_',' '), 'FontSize',10)
xlabel('Reactions', 'FontSize',xyLabelSize);
ylabel('Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 1 (by int)'], 'FontSize',titleSize)
c = colorbar; c.Label.String = 'Correlation';
c.Label.FontSize = xyLabelSize; c.FontSize = axesLabelSize;
saveas(gcf,[homePath saveFolder 'corr_int_2models_model1_rxns'],'png'); saveas(gcf,[homePath saveFolder 'corr_int_2models_model1_rxns'],'fig');

% Int - Model 2
c_int2m2_sparseCon = corr(model2{2}.int(:,sparseCon_idx_2m_noUnconstrained));
c_int2m2_rxns = corr(model2{2}.int(:,sparseCon_idx_2m_noUnconstrained)');

figure(69)
imagesc(c_int2m2_sparseCon)
set(gca, 'XTick',1:numel(model2{2}.sparse_con(sparseCon_idx_2m_noUnconstrained)), 'XTickLabel',cellstr(int2str(model2{2}.sparse_con(sparseCon_idx_2m_noUnconstrained)')), ...
    'YTick',1:numel(model2{2}.sparse_con(sparseCon_idx_2m_noUnconstrained)), 'YTickLabel',cellstr(int2str(model2{2}.sparse_con(sparseCon_idx_2m_noUnconstrained)')))
xlabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 2 (by int)'], 'FontSize',titleSize)
c = colorbar; c.Label.String = 'Correlation';
c.Label.FontSize = xyLabelSize; c.FontSize = axesLabelSize;
saveas(gcf,[homePath saveFolder 'corr_int_2models_model2_sparseCon'],'png'); saveas(gcf,[homePath saveFolder 'corr_int_2models_model2_sparseCon'],'fig');

figure(70)
imagesc(c_int2m2_rxns);
set(gca, 'XTick',1:numel(model2{2}.intl_idx), 'XTickLabel',strrep(model2{2}.rxns(model2{2}.intl_idx),'_',' '), 'XTickLabelRotation',45, ...
    'YTick',1:numel(model2{2}.intl_idx), 'YTickLabel',strrep(model2{2}.rxns(model2{2}.intl_idx),'_',' '), 'FontSize',10)
xlabel('Reactions', 'FontSize',xyLabelSize);
ylabel('Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 2 (by int)'], 'FontSize',titleSize)
c = colorbar; c.Label.String = 'Correlation';
c.Label.FontSize = xyLabelSize; c.FontSize = axesLabelSize;
saveas(gcf,[homePath saveFolder 'corr_int_2models_model2_rxns'],'png'); saveas(gcf,[homePath saveFolder 'corr_int_2models_model2_rxns'],'fig');

% Flux
c_flux2_sparseCon = corr(model2{1}.flux(model2{1}.intl_idx,sparseCon_idx_2m_noUnconstrained),model2{2}.flux(model2{2}.intl_idx,sparseCon_idx_2m_noUnconstrained));
c_flux2_rxns = corr(model2{1}.flux(model2{1}.intl_idx,sparseCon_idx_2m_noUnconstrained)',model2{2}.flux(model2{2}.intl_idx,sparseCon_idx_2m_noUnconstrained)');

figure(71)
imagesc(c_flux2_sparseCon)
set(gca, 'XTick',1:numel(model2{2}.sparse_con(sparseCon_idx_2m_noUnconstrained)), 'XTickLabel',cellstr(int2str(model2{2}.sparse_con(sparseCon_idx_2m_noUnconstrained)')), ...
    'YTick',1:numel(model2{1}.sparse_con(sparseCon_idx_2m_noUnconstrained)), 'YTickLabel',cellstr(int2str(model2{1}.sparse_con(sparseCon_idx_2m_noUnconstrained)')))
xlabel('Constraint Bound for the Number of Internal Reactions (Model 2)', 'FontSize',xyLabelSize)
ylabel('Constraint Bound for the Number of Internal Reactions (Model 1)', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models (by flux)'], 'FontSize',titleSize)
c = colorbar; c.Label.String = 'Correlation';
c.Label.FontSize = xyLabelSize; c.FontSize = axesLabelSize;
saveas(gcf,[homePath saveFolder 'corr_flux_2models_sparseCon'],'png'); saveas(gcf,[homePath saveFolder 'corr_flux_2models_sparseCon'],'fig');

figure(72)
imagesc(c_flux2_rxns);
set(gca, 'XTick',1:numel(model2{1}.intl_idx), 'XTickLabel',strrep(model2{1}.rxns(model2{1}.intl_idx),'_',' '), 'XTickLabelRotation',45, ...
    'YTick',1:numel(model2{2}.intl_idx), 'YTickLabel',strrep(model2{2}.rxns(model2{2}.intl_idx),'_',' '), 'FontSize',10)
xlabel('Reactions (Model 2)', 'FontSize',xyLabelSize);
ylabel('Reactions (Model 1)', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models (by flux)'], 'FontSize',titleSize)
c = colorbar; c.Label.String = 'Correlation';
c.Label.FontSize = xyLabelSize; c.FontSize = axesLabelSize;
saveas(gcf,[homePath saveFolder 'corr_flux_2models_rxns'],'png'); saveas(gcf,[homePath saveFolder 'corr_flux_2models_rxns'],'fig');

% Flux - Model 1
c_flux2m1_sparseCon = corr(model2{1}.flux(model2{1}.intl_idx,sparseCon_idx_2m_noUnconstrained));
c_flux2m1_rxns = corr(model2{1}.flux(model2{1}.intl_idx,sparseCon_idx_2m_noUnconstrained)');

figure(73)
imagesc(c_flux2m1_sparseCon)
set(gca, 'XTick',1:numel(model2{1}.sparse_con(sparseCon_idx_2m_noUnconstrained)), 'XTickLabel',cellstr(int2str(model2{1}.sparse_con(sparseCon_idx_2m_noUnconstrained)')), ...
    'YTick',1:numel(model2{1}.sparse_con(sparseCon_idx_2m_noUnconstrained)), 'YTickLabel',cellstr(int2str(model2{1}.sparse_con(sparseCon_idx_2m_noUnconstrained)')))
xlabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 1 (by flux)'], 'FontSize',titleSize)
c = colorbar; c.Label.String = 'Correlation';
c.Label.FontSize = xyLabelSize; c.FontSize = axesLabelSize;
saveas(gcf,[homePath saveFolder 'corr_flux_2models_model1_sparseCon'],'png'); saveas(gcf,[homePath saveFolder 'corr_flux_2models_model1_sparseCon'],'fig');

figure(74)
imagesc(c_flux2m1_rxns);
set(gca, 'XTick',1:numel(model2{1}.intl_idx), 'XTickLabel',strrep(model2{1}.rxns(model2{1}.intl_idx),'_',' '), 'XTickLabelRotation',45, ...
    'YTick',1:numel(model2{1}.intl_idx), 'YTickLabel',strrep(model2{1}.rxns(model2{1}.intl_idx),'_',' '), 'FontSize',10)
xlabel('Reactions', 'FontSize',xyLabelSize);
ylabel('Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 1 (by flux)'], 'FontSize',titleSize)
c = colorbar; c.Label.String = 'Correlation';
c.Label.FontSize = xyLabelSize; c.FontSize = axesLabelSize;
saveas(gcf,[homePath saveFolder 'corr_flux_2models_model1_rxns'],'png'); saveas(gcf,[homePath saveFolder 'corr_flux_2models_model1_rxns'],'fig');

% Flux - Model 2
c_flux2m2_sparseCon = corr(model2{1}.flux(model2{2}.intl_idx,sparseCon_idx_2m_noUnconstrained));
c_flux2m2_rxns = corr(model2{2}.flux(model2{2}.intl_idx,sparseCon_idx_2m_noUnconstrained)');

figure(75)
imagesc(c_flux2m2_sparseCon)
set(gca, 'XTick',1:numel(model2{2}.sparse_con(sparseCon_idx_2m_noUnconstrained)), 'XTickLabel',cellstr(int2str(model2{2}.sparse_con(sparseCon_idx_2m_noUnconstrained)')), ...
    'YTick',1:numel(model2{2}.sparse_con(sparseCon_idx_2m_noUnconstrained)), 'YTickLabel',cellstr(int2str(model2{2}.sparse_con(sparseCon_idx_2m_noUnconstrained)')))
xlabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
ylabel('Constraint Bound for the Number of Internal Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 2 (by flux)'], 'FontSize',titleSize)
c = colorbar; c.Label.String = 'Correlation';
c.Label.FontSize = xyLabelSize; c.FontSize = axesLabelSize;
saveas(gcf,[homePath saveFolder 'corr_flux_2models_model2_sparseCon'],'png'); saveas(gcf,[homePath saveFolder 'corr_flux_2models_model2_sparseCon'],'fig');

figure(76)
imagesc(c_flux2m2_rxns);
set(gca, 'XTick',1:numel(model2{2}.intl_idx), 'XTickLabel',strrep(model2{2}.rxns(model2{2}.intl_idx),'_',' '), 'XTickLabelRotation',45, ...
    'YTick',1:numel(model2{2}.intl_idx), 'YTickLabel',strrep(model2{2}.rxns(model2{2}.intl_idx),'_',' '), 'FontSize',10)
xlabel('Reactions', 'FontSize',xyLabelSize);
ylabel('Reactions', 'FontSize',xyLabelSize)
title(['2 ' model_title ' Models: Model 2 (by flux)'], 'FontSize',titleSize)
c = colorbar; c.Label.String = 'Correlation';
c.Label.FontSize = xyLabelSize; c.FontSize = axesLabelSize;
saveas(gcf,[homePath saveFolder 'corr_flux_2models_model2_rxns'],'png'); saveas(gcf,[homePath saveFolder 'corr_flux_2models_model2_rxns'],'fig');

end

