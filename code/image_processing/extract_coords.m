function [x0,x1,x11,x22,y0,y1,y11,y22]=extract_coords(counts,x,level)
counts1 = counts(5:end);
x0 = x(5:end);

t = find(counts1 == max(counts1));

x11 = x0(t);
y11 = counts1(t);

temp = counts1(counts1>0);
t1 = find(counts1 == min(temp(t:end)));

x22 = x0(t1);
y22 = counts1(t1);

temp = find(round(x,2)==round(level,2), 1 );
x0 = x(temp);
y0 = counts(temp);

slope = (y22 - y11) / (x22 - x11);
intercept = y11 - slope * x11;

% Calculate the slope of the perpendicular line
perpendicular_slope = -1 / slope;

% Calculate the equation of the perpendicular line using the point-slope form
perpendicular_intercept = y0 - perpendicular_slope * x0;

% Find the intersection point of the two lines
x1 = (intercept - perpendicular_intercept) / (perpendicular_slope - slope);
y1 = perpendicular_slope * x1 + perpendicular_intercept;
