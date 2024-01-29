function width = fwhm(x,y)
% function width = fwhm(x,y)
%
% Full-Width at Half-Maximum (FWHM) of the waveform y(x)
x = x(6:end);
y = y(6:end);

%Find the half max value.
halfMax = max(y) / 2;

% Find where the data first drops below half the max.
index1 = find(y >= halfMax, 1, 'first');
% Find where the data last rises above half the max.
index2 = find(y >= halfMax, 1, 'last');

fwhmx = x(index2) - x(index1);
width = index2-index1 + 1; % FWHM in indexes

end
