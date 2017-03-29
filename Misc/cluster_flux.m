function [P1,P2,axes_handle] = cluster_flux(flux_data,xnames,ynames,fig_handle)
%CLUSTER_FLUX Bicluster flux data and create heat map
%
%REQUIRED INPUTS
% flux_data: reaction fluxes, matrix [rxns x sparse_con]
% xnames: names of the x-axis variables, cell string
% ynames: names of the y-axis variables, cell string
% fig_handle: figure handles for:
%   fig_handle(1): dendrogram
%   fig_handle(2): heat map
%
%OUTPUT
% P1,P2: permutation of node labels
%   1: y-axis
%   2: x-axis
% axes_handle: axes handles for:
%   axes_handle(1): dendrogram subplot 1
%   axes_handle(2): dendrogram subplot 2
%   axes_handle(3): heat map
%
% Meghan Thommes 02/28/2017
% Based on biclustering code from David, which he got from Daniel

%%

% Cluster by Columns
D1 = pdist(flux_data); % get pairwise distance between pairs of objects
Z1 = linkage(D1); % get tree of hierarchical clusters (cluster closest distances)
% Cluster by Rows
D2 = pdist(flux_data'); % get pairwise distance between pairs of objects
Z2 = linkage(D2); % get tree of hierarchical clusters (cluster closest distances)

% Dendrograms
figure(fig_handle(1))
axes_handle(1) = subplot(2,1,1);
[~, ~, P1] = dendrogram(Z1,0, 'Labels',ynames); % create a dendrogram plot with 0 leaf nodes, get the order of node labels (top to bottom)
set(axes_handle(1), 'XTickLabelRotation',45);
ylabel('Euclidean Distance'); box on;
axes_handle(2) = subplot(2,1,2);
[~, ~, P2] = dendrogram(Z2,0, 'Labels',xnames); % create a dendrogram plot with 0 leaf nodes, get the order of node labels (top to bottom)
ylabel('Euclidean Distance'); box on;
set(axes_handle(2), 'XTickLabelRotation',45);

% Heat Map
figure(fig_handle(2))
imagesc(flux_data(P1,P2));
axes_handle(3) = gca;
set(axes_handle(3), 'XTick',1:numel(P2), 'XTickLabel',xnames(P2))
set(axes_handle(3), 'YTick',1:numel(P1), 'YTickLabel',ynames(P1))


end

