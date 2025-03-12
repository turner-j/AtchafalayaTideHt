% +------------------------------------------------------+
% |    Convert degrees to cardinal compass directions    | 
% |              with MATLAB Implementation              | 
% |                                                      |
% | Author: Ph.D. Eng. Hristo Zhivomirov        06/25/24 | 
% +------------------------------------------------------+
% 
% function: cardinal = deg2cardinal(degree)
%
% Input:
% degree - compass directions in degrees (e.g., [0 45 67.5]).
% 
% Output: 
% cardinal - converted cardinal compass directions (e.g., {'N', 'NE', 'ENE'}).
%
% Note: The function uses 16 cardinal compass directions:
%       'N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 
%       'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW'.
% If one desires only the 8 principal compass directions (N, NE, E, etc.),
% one should remove the three-letter directions /lines 34 and 35/ without
% changing any other parts of the function. Similarly, if one wants all 32
% compass directions, one should insert the quarter winds (NbE, NEbN, NEbE,
% etc.) without changing any other part of the function.

function cardinal = deg2cardinal(degree)

% input validation
validateattributes(degree, {'single', 'double'}, ...
                           {'vector', 'real', 'nonnan', 'nonempty', ...
                            '>=', 0, '<=', 360}, ...
                            '', 'degree', 1)

% convert the compass direction from deg to cardinal
CardDir = {'N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', ...
           'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW', ...
           'N'};
CardInd = discretize(degree, [0, 11.25:22.5:348.75, 360]);
cardinal = CardDir(CardInd);

end