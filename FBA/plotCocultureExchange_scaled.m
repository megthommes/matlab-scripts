function [fig_handle,axes_handle,met_handle] = plotCocultureExchange_scaled(model1_flux,model2_flux,met_conc,dt,arrow_scale,h,c,fluxType)
%PLOTCOCULTUREEXCHANGE Plots the flux of model 2 versus the flux of model 1
%with arrows pointing in the direction of flux change from time point t to
%time point t+dt. Each time point is scaled by the metabolite concentration
%   [fig_handle,axes_handle,mets_handle] = PLOTCOCULTUREEXCHANGE(model1_flux,model2_flux)
%   [fig_handle,axes_handle,mets_handle] = PLOTCOCULTUREEXCHANGE(model1_flux,model2_flux,dt,arrow_scale,h,c,fluxType)
%
% REQUIRED INPUTS
% model1_flux, model2_flux: matrix of extracellular metabolite flux over time [time x 1]
% met_conc: metabolite amount over time [time x 1]
%
% OPTIONAL INPUTS
% dt: time step as the number of indices (default=1)
% arrow_scale: length of each arrow (default=1)
% h: figure handle
% c: color matrix [1 x 3] (default=red)
% fluxType: this overrides c if included
%   0: not used by model 1 or model 2
%   1: produced by model 1 and model 2
%   2: consumed by model 1 and model 2
%   3: produced by model 1 and not used by model 2
%   4: not used by model 1 and produced by model 2
%   5: consumed by model 1 and not used by model 2
%   6: not used by model 1 and consumed by model 2
%   7: consumed by model 1 and produced by model 2
%   8: produced by model 1 and consumed by model 2
%
% OUTPUTS
% fig_handle: figure handle
% axes_handle: axes handle
% met_handle: handle for the metabolite
%
% Meghan Thommes 05/04/2017 - Based off of plotCocultureExchange

%% Check Input Variables

% Make sure model1_flux and model2_flux are matrices
if ~ismatrix(model1_flux) || ~ismatrix(model2_flux)
    error('plotCocultureExchange:incorrectInput','Error! Flux matrices must be matrices.')
end
% Make sure flux matrices are the same size
if size(model1_flux) ~= size(model2_flux)
    error('plotCocultureExchange:incorrectInput','Error! Flux matrices must be the same size.')
end
% Check for optional inputs
if nargin == 3
    dt = 1;
    arrow_scale = 1;
elseif nargin == 4
    arrow_scale = 1;
elseif nargin < 3
    error('plotCocultureExchange:incorrectInput','Error! Not enough input variables.')
end

%% Plot

if ~exist('h','var')
    fig_handle = figure;
else
    fig_handle = figure(h);
end

if ~exist('c','var')
    c = [1,0,0]; % red
end

cmap = [[0.90,0.90,0.90]; ... % 0:light gray (not used by either model)
    [0.00,0.30,0.10]; ... % 1:dark green (produced by both models)
    [0.25,0.55,0.10]; ... % 2:model2{1}.medium green (produced by model 1)
    [0.50,0.80,0.10]; ... % 3:light green (produced by model 2)
    [0.00,0.10,0.30]; ... % 4:dark blue (consumed by both models)
    [0.00,0.35,0.55]; ... % 5:model2{1}.medium blue (consumed by model 1)
    [0.00,0.60,0.80]; ... % 6:light blue (consumed by model 2)
    [1.00,0.80,0.10]; ... % 7:golden yellow (2 -> 1)
    [0.90,0.60,0.10]]; ... % 8:brown yellow (1 -> 2)

% x- and y- axes
x_limits = [floor(min(model1_flux(:))); ceil(max(model1_flux(:)))];
if abs(diff(x_limits)) <= 1E-3; x_limits = [x_limits(1)-0.5; x_limits(2)+0.5]; end
y_limits = [floor(min(model2_flux(:))); ceil(max(model2_flux(:)))];
if abs(diff(y_limits)) <= 1E-3; y_limits = [y_limits(1)-0.5; y_limits(2)+0.5]; end

markerSize = 100.*met_conc(1:dt:numel(model1_flux)-dt)./max(met_conc(1:dt:numel(model1_flux)-dt));
markerSize(isnan(markerSize)) = 1;
markerSize(find(markerSize == 0)) = 1;
met_handle = plot(model1_flux(1:dt:numel(model1_flux)-dt),model2_flux(1:dt:numel(model1_flux)-dt),'-','LineWidth',3, 'Color',c); hold on
scatter(model1_flux(1:dt:numel(model1_flux)-dt),model2_flux(1:dt:numel(model1_flux)-dt),...
    markerSize,c,'filled');
for tt = 1:dt:numel(model1_flux)-dt % each time point
    % Line between time time points 1 and 2
    x = model1_flux([tt, tt+dt]); % x-axis
    y = model2_flux([tt, tt+dt]); % y-axis
    
    % Arrow between time points 1 and 2
    p1 = [x(1),y(1)]; % first point
    mid_p = p1 + [diff(x),diff(y)]./2; % midway point
    v = [diff(x), diff(y)]; % vector
    L = sqrt(sum(v.^2)); % vector length
    if L >= arrow_scale % only plot if vector is long enough
        v_norm = arrow_scale.*v./L; % normalized vector
        quiver(mid_p(1),mid_p(2), v_norm(1),v_norm(2), 'LineWidth',3, 'Color',c, 'MaxHeadSize',1)
    end
end
hold off
xlabel('Model 1 Exchange Flux (mmol/gCDW/hr)')
ylabel('Model 2 Exchange Flux (mmol/gCDW/hr)')
axis tight

axes_handle = gca;

end

