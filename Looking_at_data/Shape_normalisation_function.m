function [outputArg4] = Shape_normalisation_function(x)
%Normalises shape of NETseq matrix of genes by subtracting the mean and 
% dividing through by the standard deviation. 
MEAN = mean(x , 2) ;
% subtract list of means from entire matrix
MINUS_MEAN = bsxfun(@minus, x, MEAN)   ;

% Find standard dev
STD = std(x, 0, 2);

% Divide through by STD
outputArg4 = bsxfun(@rdivide, MINUS_MEAN, STD);

end

