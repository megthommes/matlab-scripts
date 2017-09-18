function [fluxType] = classifyExchangeFlux(model1_flux,model2_flux,tol)
%ClassifyExchange Classify exchange flux between two models.
%   [fluxType] = ClassifyExchange(model1_flux,model2_flux)
%   [fluxType] = ClassifyExchange(model1_flux,model2_flux,tol)
%
%By default, classifyExchangeFlux assumes that exchange reactions are written
%as export reactions: 'A <==>'; therefore, uptaking a metabolite is a
%negative flux and the stoichiometric coefficient of metabolite A is -1.
%
% REQUIRED INPUTS
% model1_flux, model2_flux: matrix of extracellular metabolite flux over time
%
% OPTIONAL INPUT
% tol: tolerance
%   if flux is less than -tolerance, it is negative
%   if flux is greater than +tolerance, it is positive
%      default = 1E-6
%
% OUTPUTS
% fluxType: Number indicates what type of flux occurs (same dimensions as
%  model1_flux and model2_flux):
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
% Meghan Thommes 08/29/2017

%% Check Input Variables

% Make sure model1_flux and model2_flux are matrices
if ~ismatrix(model1_flux) || ~ismatrix(model2_flux)
    error('classifyExchangeFlux:IncorrectInput', ...
        'Error! Flux matrices must be matrices.')
end
% Make sure flux matrices are the same size
if size(model1_flux) ~= size(model2_flux)
    error('classifyExchangeFlux:IncorrectInput', ...
        'Error! Flux matrices must be the same size.')
end
% Check for optional inputs
if nargin < 2
    error('classifyExchangeFlux:IncorrectInput', ...
        'Error! Not enough input variables.')
end

if ~exist('tol','var')
    tol = 1E-6;
end

%% Classify

% Determine if metabolites are produced, consumed, or not used
produced1 = find(model1_flux > tol); % produced by model 1
produced2 = find(model2_flux > tol); % produced by model 2
consumed1 = find(model1_flux < -tol); % consumed by model 1
consumed2 = find(model2_flux < -tol); % consumed by model 2
notused1 = find(model1_flux > -tol & model1_flux < tol); % not used by model 1
notused2 = find(model2_flux > -tol & model2_flux < tol); % not used by model 2

% Find the indices
idx_notused1_notused2 = intersect(notused1,notused2); % neither produced nor consumed
idx_pro1_pro2 = intersect(produced1,produced2); % produced by model 1 & model 2
idx_con1_con2 = intersect(consumed1,consumed2); % consumed by model 1 & model 2
idx_pro1_notused2 = intersect(produced1,notused2); % produced by model 1 & not used by model 2
idx_notused1_pro2 = intersect(produced2,notused1); % produced by model 2 & not used by model 1
idx_con1_notused2 = intersect(consumed1,notused2); % consumed by model 1 & not used by model 2
idx_notused1_con2 = intersect(consumed2,notused1); % consumed by model 2 & not used by model 1
idx_con1_pro2 = intersect(produced2,consumed1); % produced by model 2 & consumed by model 1
idx_pro1_con2 = intersect(produced1,consumed2); % produced by model 1 & consumed by model 2

% Classify exchange flux
fluxType = zeros(size(model1_flux));
fluxType(idx_notused1_notused2) = 0;
fluxType(idx_pro1_pro2) = 1;
fluxType(idx_pro1_notused2) = 2;
fluxType(idx_notused1_pro2) = 3;
fluxType(idx_con1_con2) = 4;
fluxType(idx_con1_notused2) = 5;
fluxType(idx_notused1_con2) = 6;
fluxType(idx_con1_pro2) = 7;
fluxType(idx_pro1_con2) = 8;

end