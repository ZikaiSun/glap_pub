% function I = histomatch(I, J, mask);
%     h1 = sort(J(:));
%     [h2, j] = sort(I(:));
%     I(j) = h1;
% end

function I = histomatch(I, J, mask)
    % mask = logical(ones(size(I)));
    % Ensure mask is of logical type
    mask = logical(mask);

    % Extract and sort values from J within the mask-specified region
    J_masked = J(mask);
    h1 = sort(J_masked(:));

    % Extract values from I within the mask-specified region and sort them
    I_masked = I(mask);
    [h2, j] = sort(I_masked(:));

    % Create an array of the same size as h2 for storing the mapped values
    mapped_values = zeros(size(h2));

    % Ensure the array for the mapped values has the same number of elements as the original
    % If there are fewer elements in J than I within the region, handle appropriately
    min_length = min(length(h1), length(h2));
    h1 = h1(1:min_length);
    h2 = h2(1:min_length);
    j = j(1:min_length);

    % Apply mapping to the indices sorted by h2
    mapped_values(1:min_length) = h1;

    % Place the mapped values back into their corresponding locations in I
    I_temp = I_masked;
    I_temp(j) = mapped_values;

    % Insert the processed region back into the original image
    I(mask) = I_temp;
end

% function I = histomatch(I, J, mask)
%     % Ensure mask is of logical type
%     mask = logical(mask);

%     % Extract and sort values from J within the mask-specified region
%     J_masked = J(mask);
%     h1 = sort(J_masked(:));

%     % Extract values from I within the mask-specified region and sort them
%     I_masked = I(mask);
%     [h2, j] = sort(I_masked(:));

%     % Create an array for the mapped values
%     mapped_values = zeros(size(h2));

%     % Determine the minimum length for mapping
%     min_length = min(length(h1), length(h2));
%     h1 = h1(1:min_length);
%     h2 = h2(1:min_length);
%     j = j(1:min_length);

%     % Apply mapping to the indices sorted by h2
%     mapped_values(1:min_length) = h1;

%     % Map the entire image I using the established mapping from the mask region 
%     % Create a scatter mapping between sorted I values from mask and corresponding J values
%     scatter_map = containers.Map('KeyType', 'double', 'ValueType', 'double');

%     % Fill the mapping with matched values from the masked area
%     for index = 1:min_length
%         scatter_map(h2(index)) = h1(index);
%     end

%     % Create a handle for missing entries: map to the nearest in scatter_map
%     keys = h2;
%     values = h1;
    
%     % Apply this mapping to each element in the entire image I
%     I_all_values = I(:);
%     for i = 1:length(I_all_values)
%         current_value = I_all_values(i);
        
%         % If present, use the mapped value; otherwise find the closest
%         if isKey(scatter_map, current_value)
%             I_all_values(i) = scatter_map(current_value);
%         else
%             % Find the closest h2 key
%             [~, closest_index] = min(abs(keys - current_value));
%             I_all_values(i) = values(closest_index);
%         end
%     end

%     % Reshape the result back to the original image size
%     I = reshape(I_all_values, size(I));
% end

% function I = histomatch(I, J, mask)
%     % Ensure mask is of logical type
%     mask = logical(mask);

%     % Extract and sort values from J within the mask-specified region
%     J_masked = J(mask);
%     h1 = sort(J_masked(:));

%     % Extract values from I within the mask-specified region and sort them
%     I_masked = I(mask);
%     [h2, j] = sort(I_masked(:));

%     % Create an array of the same size as h2 for storing the mapped values
%     mapped_values = zeros(size(h2));

%     % Ensure the array for the mapped values has the same number of elements as the original
%     % If there are fewer elements in J than I within the region, handle appropriately
%     min_length = min(length(h1), length(h2));
%     h1 = h1(1:min_length);
%     h2 = h2(1:min_length);

%     % Apply mapping to the indices sorted by h2
%     mapped_values(1:min_length) = h1;

%     % Create a lookup table for value mapping based on h2 and mapped_values
%     value_map = containers.Map(h2, mapped_values);

%     % Handle values outside the masked region using nearest neighbor interpolation
%     % Complete mapping including all unique values in I
%     unique_values_in_I = unique(I(:));

%     % Create extended value map covering all unique values from the original image
%     extended_value_map = zeros(size(unique_values_in_I));
%     for k = 1:length(unique_values_in_I)
%         val = unique_values_in_I(k);
%         if isKey(value_map, val)
%             extended_value_map(k) = value_map(val);
%         else
%             % Nearest neighbor lookup for unmapped values
%             [~, idx] = min(abs(h2 - val));
%             extended_value_map(k) = mapped_values(idx);
%         end
%     end

%     % Apply this extended mapping to the entire image
%     % Instead of I + 1, directly use the value map to substitute values
%     for k = 1:length(unique_values_in_I)
%         I(I == unique_values_in_I(k)) = extended_value_map(k);
%     end
% end






