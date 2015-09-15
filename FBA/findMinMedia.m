
clear all
close all
clc

yeast7_1 = load('C:\Users\mthommes\Dropbox\Segre\Epistasis\Yeast\YeastNet_Models\Version7\yeast7_1.mat');
yeast6 = load('C:\Users\mthommes\Dropbox\Segre\Epistasis\Yeast\YeastNet_Models\Version6\yeast6_06.mat');
yeast7 = load('C:\Users\mthommes\Dropbox\Segre\Epistasis\Yeast\YeastNet_Models\Version7\yeast7_00.mat');
initCobraToolbox

model = yeast6;

%% Minimal Media

clc

% start with a clean slate: unconstrained excretion, no uptake)
exchangeRxns = findExcRxns(model);
model.lb(exchangeRxns) = -1;
model.ub(exchangeRxns) = 1000;
desiredExchanges = {... % yeast6
    'r_1654'; ... % 'ammonium exchange'
    'r_1832'; ... % 'H+ exchange'
    'r_1861'; ... % 'iron(2+) exchange'
    'r_1992'; ... % 'oxygen exchange'
    'r_2005'; ... % 'phosphate exchange'
    'r_2049'; ... % 'sodium exchange'
    'r_2060'; ... % 'sulphate exchange'
    'r_2100'; ... % 'water exchange'
    };
% desiredExchanges = {... % yeast7
%     'r_1654'; ... % 'ammonium exchange'
%     'r_1832'; ... % 'H+ exchange'
%     'r_1861'; ... % 'iron(2+) exchange'
%     'r_1992'; ... % 'oxygen exchange'
%     'r_2005'; ... % 'phosphate exchange'
%     'r_2020'; ... % 'potassium exchange'
%     'r_2049'; ... % 'sodium exchange'
%     'r_2060'; ... % 'sulphate exchange'
%     'r_2100'; ... % 'water exchange'
%     };
glucoseExchange = {...
    'r_1714'; ... % D-glucose exchange'
    };

uptakeRxnIndices = findRxnIDs(model,desiredExchanges);
glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);

% modify uptake
model.lb(uptakeRxnIndices) = -1000;
% model.lb(glucoseExchangeIndex) = -10;
sol = optimizeCbModel(model)

return
%% Save COMETS Model

exRxns = find(exchangeRxns == 1);
cobra2comets_new('yeast6_06.comets',model,find(model.c==1),exRxns);

%%
clc
model = yeast7_1;
minMedia = {...
    's_0420'; ... % 'ammonium [extracellular]'
    's_0565'; ... % 'D-glucose [extracellular]'
    's_0796'; ... % 'H+ [extracellular]'
    's_0925'; ... % 'iron(2+) [extracellular]'
    's_1277'; ... % 'oxygen [extracellular]'
    's_1324'; ... % 'phosphate [extracellular]'
    's_1438'; ... % 'sodium [extracellular]'
    's_1468'; ... % 'sulphate [extracellular]'
    's_0805'; ... % 'H2O [extracellular]'
    's_1374'; ... % 'potassium [extracellular]'
    };
minMediaIndices = findMetIDs(model,minMedia);
[model.metNames(minMediaIndices) minMedia]

%%
clc
model = yeast7_1;
exchangeRxns = findExcRxns(model);
exRxns = find(exchangeRxns == 1);

[exchangeRxns,nutrients] = findExcRxns(model);
exRxns = find(exchangeRxns == 1);
temp = model.S(:,exRxns);
[exMets,temp2] = find(temp ~= 0);
[model.metNames(exMets) model.mets(exMets)]
whos exRxns