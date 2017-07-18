function [letter] = num2letter(num)
%NUM2LETTER Convert a number into a letter (eg 1 to A, 27 to AA)
%
% Only valid for numbers up to and including 208
%
%REQUIRED INPUT
% num
%
%OUTPUT
% letter
%
% Meghan Thommes 07/13/2017 -  Up to and including 208 instead of 104
% Meghan Thommes 05/16/2017

%% Check Inputs

if size(num,1) < size(num,2)
    num = num';
end

%%

if num <= 26
    letter = char(num-1+'A');
elseif num > 26 && num <= 52
    letter = ['A' char(num-27+'A')];
elseif num > 52 && num <= 78
    letter = ['B' char(num-53+'A')];
elseif num > 78 && num <= 104
    letter = ['B' char(num-79+'A')];
elseif num > 104 && num <= 130
    letter = ['C' char(num-105+'A')];
elseif num > 130 && num <= 156
    letter = ['D' char(num-131+'A')];
elseif num > 156 && num <= 182
    letter = ['E' char(num-157+'A')];
elseif num > 182 && num <= 208
    letter = ['F' char(num-183+'A')];
end

% letter = char(num-1+'A'-26*(ceil(num./26)-1));


end

