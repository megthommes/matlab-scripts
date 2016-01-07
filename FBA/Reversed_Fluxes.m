function [new_model] = Reversed_Fluxes(model)
%Reversed_Fluxes Reverse the stoichiometric coefficients so that the
%exchange reactions are -1 instead of +1
%model = Reversed_Fluxes(model)
%
%REQUIRED INPUT
% The model structure must contain the following fields:
%   model.S:     Stoichiometric matrix
%   model.c:     Objective coefficients
%   model.lb:    Lower bounds
%
%OUTPUT
% model is a new model with the updated stoichiometric matrix
%
%3/4/2015 - Meghan Thommes

disp(strrep(model.description,'_',' '))

% Exchange and Nutrient Reactions
[exchRxns,nutRxns] = findExchRxns(model,1);
[fakeExchRxns,~] = findExchRxns(model);
% exchRxns(1620) = false; % 2,3-cyclic AMP production - Yarrowia
new_model = model;
disp([num2str(sum(exchRxns)) ' exchange reactions & ' num2str(sum(nutRxns)) ' nutrients'])
new_model.S(:,exchRxns) = -1.*model.S(:,exchRxns);
new_model.lb(exchRxns) = -1.*model.ub(exchRxns);
new_model.ub(exchRxns) = -1.*model.lb(exchRxns);
new_model.S(:,fakeExchRxns) = -1.*model.S(:,fakeExchRxns);
new_model.lb(fakeExchRxns) = -1.*model.ub(fakeExchRxns);
new_model.ub(fakeExchRxns) = -1.*model.lb(fakeExchRxns);
fprintf('\n')
% Objective Value
sol = FBA(new_model);
disp(['Growth Flux is ' num2str(sol.objectiveValue)])
disp(['Status is ' sol.status])

end