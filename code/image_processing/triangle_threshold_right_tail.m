function [threshold, line_endpoints] = triangle_threshold_right_tail(histogram)
    % Find the maximum value and its location in the histogram
    [max_value, max_index] = max(histogram);
    
    % Find the location of the minimum value on the right side
    [~, min_index] = min(histogram(max_index:end));
    min_index = min_index + max_index - 1;
    
    % Calculate the slope and intercept of the line
    slope = (histogram(min_index) - max_value) / (min_index - max_index);
    intercept = max_value - slope * max_index;
    
    % Construct the line endpoints
    line_endpoints = [max_index, max_value; min_index, histogram(min_index)];
    
    % Construct the line as an array of y values
    line = slope * (1:(min_index - max_index + 1)) + intercept;
    
    % Calculate the distances between the histogram and the line
    distances = line-histogram(max_index:min_index);
    
    % Find the threshold value as the index with the maximum distance
    [~, threshold_index] = max(distances);
    
    % Calculate the threshold value
    threshold = max_index + threshold_index - 1;
end