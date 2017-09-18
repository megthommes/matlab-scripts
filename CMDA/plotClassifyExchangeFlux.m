function [fig_handle,axes_handle,cbar_handle] = plotClassifyExchangeFlux(fluxType,fig_h,axes_h,metNames)
%plotClassifyExchangeFlux Plots the type of metabolite exchange between two
%models. Uses the output from classifyExchangeFlux.
%   [fig_handle,axes_handle,cbar_handle] = plotClassifyExchangeFlux(fluxType)
%   [fig_handle,axes_handle,cbar_handle] = plotClassifyExchangeFlux(fluxType,fig_h,axes_h,metNames)
%
%REQUIRED INPUTS
% fluxType: Number indicates what type of flux occurs [metabolites x time]
%   0: not used by model 1 or model 2 (->)
%   1: produced by model 1 and model 2 (1,2->)
%   2: produced by model 1 and not used by model 2 (1->)
%   3: not used by model 1 and produced by model 2 (2->)
%   4: consumed by model 1 and model 2 (->1,2)
%   5: consumed by model 1 and not used by model 2 (->1)
%   6: not used by model 1 and consumed by model 2 (->2)
%   7: consumed by model 1 and produced by model 2 (2->1)
%   8: produced by model 1 and consumed by model 2 (1->2)
%
%OPTIONAL INPUTS
% fig_h: figure handle
% axes_h: axes handle
% metNames: names of the metabolites
%
%OUTPUTS
% fig_handle: figure handle
% axes_handle: axes handle
% cbar_handle: colorbar handle
%
% Meghan Thommes 08/29/2017

%% Check Input Variables

% Make sure have enough inputs
if nargin < 1
    error('plotClassifyExchangeFlux:incorrectInput', ...
        'Error! Not enough input variables.')
end
% Check for optional inputs
if ~exist('fig_h','var') || isempty(fig_h)
    fig_handle = figure;
else
    fig_handle = figure(fig_h);
end
if ~exist('axes_h','var') || isempty(axes_h)
    axes_handle = gca;
else
    axes_handle = axes_h;
end

%% Plot

% Load Colormap Data
load('classifyExchangeFlux_cmap.mat')

imagesc(axes_handle,fluxType)
colormap(fluxType_cmap)
caxis([0 8])
% colorbar labels
cbar_handle = colorbar;
cbar_handle.Label.String = 'Exchange Type';
cbar_handle.Ticks = 0.45:8/9:8;
cbar_handle.TickLabels = {'\rightarrow','1,2\rightarrow','1\rightarrow','2\rightarrow',...
    '\rightarrow1,2','\rightarrow1','\rightarrow2','2\rightarrow1','1\rightarrow2'};
% y-axis labels
if exist('metNames','var')
    set(axes_handle, 'YTick',1:numel(metNames), 'YTickLabel',metNames)
end

ylabel('Extracellular Metabolites')

end

