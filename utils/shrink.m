function smask = shrink(mask, t)
    % Shrink the mask, it uses the imdilate function, so the periphery must be 0
    % TODO: Set the outermost row and column (i.e., the edges) of the mask to 0
    
    % Set the border of the mask to 0
    mask(1,:) = 0;            % Set the first row to 0
    mask(:,1) = 0;            % Set the first column to 0
    mask(end,:) = 0;          % Set the last row to 0
    mask(:,end) = 0;          % Set the last column to 0
    
    % Create a structuring element
    se = strel('disk', t);    % Originally was 'se = ones(t, t);'
    
    % Perform dilation on the inverted mask and then invert it back
    smask = imdilate(1 - mask, se);
    smask = logical(1 - smask);
end