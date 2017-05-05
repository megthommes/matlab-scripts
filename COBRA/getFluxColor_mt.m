function color = getFluxColor_mt(fluxValues,fluxScale)
%GETFLUXCOLOR_MT Obtains color 
%
% color = GETFLUXCOLOR_MT(fluxValues,fluxScale)
%
%INPUT
% fluxValues:   Vector containing values
%
%OPTIONAL INPUT
% fluxScale:    Min/max flux value (default=[-1000,1000])
%
%OUTPUT
% color         n x 3 matrix containing RGB color values for each input
%               value
%
% Meghan Thommes 2/17/2017
% Richard Que 12/2009

% -Inf: [0,0,255]
% -0.01: [204,204,255]
% 0: [0.8,0.8,0.8]
% 0.01: [255,204,204]
% Inf: [255,0,0]

if nargin <2
    fluxScale = [-1000, 1000];
end

sample_flux = [fluxScale(1),-0.01,0,0.01,fluxScale(2)]';
sample_color = [[0,0,255];[204,204,255];[204,204,204];[255,204,204];[255,0,0]];
color = interp1(sample_flux,sample_color,fluxValues);

end