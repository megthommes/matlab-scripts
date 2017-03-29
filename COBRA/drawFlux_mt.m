function options = drawFlux_mt(map,model,flux,options,varargin)
%DRAWFLUX_MT overlays a flux distribution onto a reaction map
%
% options = DRAWFLUX_MT(map,model,flux,options,varargin)
%
%INPUTS
% map               map structure
% model             COBRA model structure
% flux              Flux vector to overlay
%
%OPTIONAL INPUTS
% Optional parameters can be set using either the options structure, a
% parameter name / value pair input arguments, or a combination of both.
%
% options            Structure containing optional parameters
%   fileName          Name of output file
%   lb                Lower limit to round smaller values up to.
%   ub                Upper limit to round larger values down to.
%   rxnDirMultiplier  scaling value of arrows denoting flux direction
%   fluxScale         Min/max flux value (default=[-1000,1000])
%
%  Note: see setMapOptions for additional options.
%
% varargin          optional parameter name / parameter value pairs
%
%OUTPUT
% options           Structure containing optional parameters.
%
% Meghan Thommes 2/17/2017 - Deleted colorScale, zeroFluxWidth, &
%                            zeroFluxColor. Modified flux so that it scales
%                            from minFlux to maxFlux.
%

%% Check Function Inputs

if nargin<4, options=[]; end
% Parse optional parameters
if mod(length(varargin),2)==0
    for i=1:2:length(varargin)-1
        options = setMapOptions(options,map,model,varargin{i},varargin{i+1});
    end
else
    error('Invalid number of parameters/values');
end

% If fields are empty, add default values
if ~isfield(options,'scaleType'), options.scaleType=1; end
if ~isfield(options,'lb'), lb=[]; else lb = options.lb; end
if ~isfield(options,'ub'), ub=[]; else ub = options.ub; end
if ~isfield(options,'rxnDirMultiplier'), options.rxnDirMultiplier = 2; end
if ~isfield(options,'rxnDirFlag'), rxnDirFlag = false; else rxnDirFlag = options.rxnDirFlag; end
if ~isfield(options,'fluxScale'), options.fluxScale = [-1000, 1000]; end
absFlag=false;
switch lower(options.scaleType)
    case {1, 'linear'}
        options.scaleTypeLabel='Linear;';
    case {2 ,'linear absolute'}
        flux=abs(flux);
        absFlag=true;
        options.scaleTypeLabel='Linear absolute;';
    case {3,'log10'}
        flux = log10(abs(flux));
        options.scaleTypeLabel='Log10;';
end

if ~isempty(ub)
    flux(flux>ub) = ub; % flux above upper bound is equal to upper bound
    options.overlayUB = [num2str(ub) '+'];
    fluxMax = ub;
else
    fluxMax = max(flux(~isinf(flux))); % max flux value excluding infinite
    options.overlayUB = num2str(fluxMax);
end

if ~isempty(lb)
    flux(flux<lb) = lb; % flux below lower bound is equal to lower bound
    options.overlayLB = [num2str(lb) '-'];
    fluxMin = lb;
elseif absFlag
    options.overlayLB = '0';
    fluxMin = 0;
else
    fluxMin = min(flux(~isinf(flux))); % min flux value excluding infinite
    options.overlayLB = num2str(fluxMin);
end

%%

%rxnDirectionality
if rxnDirFlag
    options.rxnDir = zeros(length(map.connectionAbb),1);
    for i = 1:length(map.connectionAbb)
        options.rxnDir(ismember(map.connectionAbb,model.rxns(flux>0))) = 1;
        options.rxnDir(ismember(map.connectionAbb,model.rxns(flux<0))) = -1;
    end
end

% Color & scale the width of the reaction arrows
options.colorScale = unique(round(getFluxColor_mt(sort(flux),options.fluxScale)),'rows','stable');
color = getFluxColor_mt(flux,options.fluxScale);
weight = getFluxWeight_mt(flux,options.fluxScale);
options.edgeWeight = zeros(size(map.connection,1),1); % width of reaction arrows
for ii = 1:numel(weight)
    idx = strcmp(map.connectionAbb,model.rxns(ii));
    options.edgeWeight(idx,1) = weight(ii); % width of the reaction arrows
    options.edgeColor(idx,:) = repmat(color(ii,:),sum(idx),1); % color of reaction arrows    
    options.edgeArrowColor(idx,:) = repmat(color(ii,:),sum(idx),1); % color of reaction arrows
end


options.lb = fluxMin;
options.ub = fluxMax;
options.overlayType = 'Flux';
drawCbMap_mt(map,options);






end