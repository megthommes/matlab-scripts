function [fig_handle,axes_handle,mets_handle] = plotCocultureExchange(model1_flux,model2_flux,dt,arrow_scale,met_idx,h)
%PLOTCOCULTUREEXCHANGE Plots the flux of model 2 versus the flux of model 1
%with arrows pointing in the direction of flux change from time point t to
%time point t+dt
%   [fig_handle,axes_handle,mets_handle] = PLOTCOCULTUREEXCHANGE(model1_flux,model2_flux)
%   [fig_handle,axes_handle,mets_handle] = PLOTCOCULTUREEXCHANGE(model1_flux,model2_flux,dt,arrow_scale,met_idx,h)
%
% REQUIRED INPUTS
% model1_flux, model2_flux: matrix of extracellular metabolite flux over time [time x metabolites]
%
% OPTIONAL INPUTS
% dt: time step as the number of indices. (default=1)
% arrow_scale: length of each arrow. (default=1)
% met_idx: index of metabolites to plot. (default=all)
% h: figure handle
%
% OUTPUTS
% fig_handle: figure handle
% axes_handle: axes handle
% mets_handle: handle for each metabolite plotted (size of met_idx)
%
% Meghan Thommes 02/02/2017

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
if nargin == 2
    dt = 1;
    arrow_scale = 1;
    met_idx = 1:size(model1_flux,2);
elseif nargin == 3
    arrow_scale = 1;
    met_idx = 1:size(model1_flux,2);
elseif nargin == 4
    met_idx = 1:size(model1_flux,2);
elseif nargin < 2
    error('plotCocultureExchange:incorrectInput','Error! Not enough input variables.')
end

%% Plot

if ~exist('h','var')
    fig_handle = figure;
else
    fig_handle = figure(h);
end

c = distinguishable_colors(numel(met_idx));
% mets_handle = 

for mm = 1:numel(met_idx) % each metabolite
    for tt = 1:dt:size(model1_flux,1)-dt % each time point
        % Line between time time points 1 and 2
        x = model1_flux([tt, tt+dt], met_idx(mm)); % x-axis
        y = model2_flux([tt, tt+dt], met_idx(mm)); % y-axis
        mets_handle(mm) = plot(x,y,'o-','LineWidth',3, 'Color',c(mm,:)); hold on
        
        % Arrow between time points 1 and 2
        p1 = [x(1),y(1)]; % first point
        mid_p = p1 + [diff(x),diff(y)]./2; % midway point
        v = [diff(x), diff(y)]; % vector
        L = sqrt(sum(v.^2)); % vector length
        if L >= arrow_scale % only plot if vector is long enough
            v_norm = arrow_scale.*v./L; % normalized vector
            quiver(mid_p(1),mid_p(2), v_norm(1),v_norm(2), 'LineWidth',3, 'Color',c(mm,:), 'MaxHeadSize',1)
        end
    end
end
hold off
xlabel('Model 1 Exchange Flux')
ylabel('Model 2 Exchange Flux')

axes_handle = gca;

end

