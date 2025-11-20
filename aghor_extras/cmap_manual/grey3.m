function J = grey3(m)
% grey colormap with only 3 values, 0, 0.5, 1
% set the "breakpoints" for the color curve:
lowValue = 0;
midValue = 128;
highValue = 255;

% pick "any" three colors to correspond to the breakpoints:
lowColor = [255 0 0];
midColor = [40 40 40];
highColor = [0 255 255];

% create the colormap:
J = interp1( [lowValue midValue highValue], ...
  [lowColor; midColor; highColor]/255, ...
  linspace(lowValue, highValue, 256));
return J

