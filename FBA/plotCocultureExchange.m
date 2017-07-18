function [fig_handle,axes_handle,met_handle] = plotCocultureExchange(model1_flux,model2_flux,dt,fig_h,axes_h,markerSize,arrowHead,lineWidth,lineColor)
%PLOTCOCULTUREEXCHANGE Plots the flux of model 2 versus the flux of model 1
%with arrows pointing in the direction of flux change from time point t to
%time point t+dt
%   [fig_handle,axes_handle,mets_handle] = PLOTCOCULTUREEXCHANGE(model1_flux,model2_flux)
%   [fig_handle,axes_handle,mets_handle] = PLOTCOCULTUREEXCHANGE(model1_flux,model2_flux,dt,fig_h,axes_h,markerSize,arrowHead,lineWidth,lineColor)
%
% REQUIRED INPUTS
% model1_flux, model2_flux: matrix of extracellular metabolite flux over time [time x 1]
%
% OPTIONAL INPUTS
% dt: time step as the number of indices (default=1)
% fig_h: figure handle
% axes_h: axes handle
% markerSize: maximum markerSize (default=250)
% arrowHead: length and width of each arrow (default=5)
% lineWidth: width of the line (default=3)
% lineColor: color of each line (default=black)
%
% OUTPUTS
% fig_handle: figure handle
% axes_handle: axes handle
% met_handle: handle for each metabolite plotted (size of met_idx)
%
% Meghan Thommes 05/05/2017 - Added x- and y-axis lines, use annotations
%                             instead of quiver to draw arrows
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
% Make sure have enough inputs
if nargin < 2
    error('plotCocultureExchange:incorrectInput','Error! Not enough input variables.')
end
% Check for optional inputs
if ~exist('dt','var')
    dt = 1;
end
if ~exist('fig_h','var')
    fig_handle = figure;
else
    fig_handle = figure(fig_h);
end
if ~exist('axes_h','var')
    axes_handle = gca;
else
    axes_handle = axes_h;
end
if ~exist('markerSize','var')
    markerSize = 250;
end
if ~exist('arrowHead','var')
    arrowHead = 1;
end
if ~exist('lineWidth','var')
    lineWidth = 3;
end
if ~exist('lineColor','var')
    lineColor = [0,0,0];
end

%% Plot

% x-axis
x_limits = [-max(abs(([model1_flux(:); model1_flux(:)]))), max(abs(([model1_flux(:); model1_flux(:)])))];
if abs(diff(x_limits)) <= 1E-3; x_limits = [x_limits(1)-1; x_limits(2)+1]; end
plot(axes_handle,x_limits(1):diff(x_limits)/2:x_limits(2),zeros(size(x_limits(1):diff(x_limits)/2:x_limits(2))),'k--', 'LineWidth',lineWidth); hold on

% y-axis
y_limits = [-max(abs(([model2_flux(:); model2_flux(:)]))), max(abs(([model2_flux(:); model2_flux(:)])))];
if abs(diff(y_limits)) <= 1E-3; y_limits = [y_limits(1)-1; y_limits(2)+1]; end
plot(axes_handle,zeros(size(y_limits(1):diff(y_limits)/2:y_limits(2))),y_limits(1):diff(y_limits)/2:y_limits(2),'k--', 'LineWidth',lineWidth);

% line
met_handle = plot(axes_handle,model1_flux(1:dt:numel(model1_flux)),model2_flux(1:dt:numel(model1_flux)),'-','LineWidth',lineWidth, 'Color',lineColor);

% markers for relative concentration
scatter(axes_handle,model1_flux(1:dt:numel(model1_flux)),model2_flux(1:dt:numel(model1_flux)),markerSize,lineColor,'filled');

% arrows
arrowLength = 0.01;
for tt = 1:dt:numel(model1_flux)-dt % each time point
    x = model1_flux([tt, tt+dt]); % x-axis
    y = model2_flux([tt, tt+dt]); % y-axis
    
    % Arrow between time points 1 and 2
    p1 = [x(1),y(1)]; % first point
    mid_p = p1 + [diff(x),diff(y)]./2; % midway point
    v = [diff(x), diff(y)]; v(abs(v)<1E-6) = 0; % vector
    L = sqrt(sum(v.^2)); % vector length
    v_norm = arrowLength.*v./L; v_norm(isnan(v_norm)) = 0; % normalized vector
    if L > 0.1*max([x_limits(2),y_limits(2)])
        ah = annotation(fig_h,'arrow', 'HeadStyle','vback2', 'HeadLength',arrowHead, 'HeadWidth',arrowHead);
        set(ah, 'Parent',gca, 'Position',[mid_p(1),mid_p(2),v_norm(1),v_norm(2)], 'Color','k')
    end
end
hold off
xlabel('Model 1 Exchange Flux')
ylabel('Model 2 Exchange Flux')
axis tight

end

