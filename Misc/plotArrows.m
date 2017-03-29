function [fig_handle,axes_handle,scatter_handle,line_handle] = plotArrows(xdata,ydata,arrow_scale,h,c,lineWidth,lineStyle,markerStyle,markerSize)
%PLOTARROWS Plots with arrows
%   [fig_handle,axes_handle,scatter_handle,line_handle] = PLOTARROWS(xdata,ydata,arrow_scale)
%   [fig_handle,axes_handle,scatter_handle,line_handle] = PLOTARROWS(xdata,ydata,arrow_scale,h,c,lineWidth,lineStyle,markerStyle,markerSize)
%
% REQUIRED INPUTS
% xdata, ydata: matrix of data to plot
%
% OPTIONAL INPUTS
% arrow_scale: length of each arrow (default=1)
% h: figure handle
% c: color matrix [num_data x 3]
% lineWidth: line width (default=3)
% lineStyle: line style (default='-')
% markerStyle: marker symbol (default='o')
% markerSize: marker size (default=50)
%
% OUTPUTS
% fig_handle: figure handle
% axes_handle: axes handle
% scatter_handle: scatter plot handle
% line_handle: plot handle
%
% Meghan Thommes 03/20/2017 - Centered the arrows
% Meghan Thommes 03/17/2017
% Based on plotCocultureExchange

%% Check Input Variables

% Make sure model1_flux and model2_flux are matrices
if ~ismatrix(xdata) || ~ismatrix(ydata)
    error('plotArrows:incorrectInput','Error! xdata and ydata must be matrices.')
end
% Make sure flux matrices are the same size
if size(xdata) ~= size(ydata)
    error('plotArrows:incorrectInput','Error!  matrices must be the same size.')
end
% Check for optional inputs
if nargin < 2
    error('plotArrows:incorrectInput','Error! Not enough input variables.')
end

if ~exist('arrow_scale','var') || isempty(arrow_scale)
    arrow_scale = 1;
end

if ~exist('h','var') || isempty(h)
    fig_handle = figure;
else
    fig_handle = figure(h);
end

if ~exist('c','var') || isempty(c)
    c = distinguishable_colors(numel(met_idx));
end

if ~exist('lineWidth','var') || isempty(lineWidth)
    lineWidth = 3;
end

if ~exist('lineStyle','var') || isempty(lineWidth)
    lineStyle = '-';
end

if ~exist('markerStyle','var') || isempty(markerStyle)
    markerStyle = 'o';
end

if ~exist('markerSize','var') || isempty(markerSize)
    markerSize = 50;
end


%% Plot

for yy = 1:size(xdata,2)
    line_handle(yy) = plot(xdata,ydata, 'LineStyle',lineStyle, 'LineWidth',lineWidth, 'Color',c(yy,:));
    scatter_handle(yy) = scatter(xdata(:,yy),ydata(:,yy),markerSize,c(yy,:),markerStyle,'filled');
    
    for xx = 1:size(xdata,1)-1
        % Line between time time points 1 and 2
        x = xdata([xx, xx+1], yy); % x-axis
        y = ydata([xx, xx+1], yy); % y-axis
                        
        % Arrow between time points 1 and 2
        p1 = [x(1),y(1)]; % first point
        v = [diff(x), diff(y)]; % vector
        mid_p = p1 + v./2; % midway point
        L = sqrt(sum(v.^2)); % vector length
        if L >= 2*arrow_scale % only plot if vector is long enough
            v_norm = arrow_scale.*v./L; % normalized vector
            % Center the arrows
            arrow_p = mid_p - v_norm./2;
            quiver(arrow_p(1),arrow_p(2), v_norm(1),v_norm(2), 'LineWidth',lineWidth, 'Color',c(yy,:), 'MaxHeadSize',1)
        end
    end
end
hold off

axes_handle = gca;

end

