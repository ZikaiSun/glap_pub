function image = readImage(imagePath, resizeDimension)
    % turn color/blackwhite image into [0, 255] double 2d array.
    % If resizeDimension is specified, resize the shorter side of the image to that value.
    %
    % Arguments:
    % - imagePath: Path to the image file
    % - resizeDimension: Desired size of the shorter side of the image (optional)
    
    if nargin < 2
        resizeDimension = -1;
    end

    % Construct the full file path
    fullImagePath = fullfile(pwd, imagePath);

    % Check if the file exists
    if ~exist(fullImagePath, 'file')
        error('File does not exist: %s', fullImagePath);
    end

    % If the file exists, read it
    rawImage = imread(fullImagePath);

    % Resize the image if resizeDimension is specified
    if resizeDimension > 0
        rawImage = resizeImage(rawImage, resizeDimension);
    end

    % Convert the image to grayscale and double
    image = double(rgb2grey1(rawImage));
end

function gray = rgb2grey1(rgb)
    if ndims(rgb) == 3
        gray = 0.2989 * rgb(:,:,1) + 0.5870 * rgb(:,:,2) + 0.1140 * rgb(:,:,3);
    else
        gray = rgb;
    end
end

function resizedImage = resizeImage(image, resizeDimension)
    [rows, cols, ~] = size(image);

    % Determine the scaling factor to resize the shorter side to resizeDimension
    if rows < cols
        scaleFactor = resizeDimension / rows;
    else
        scaleFactor = resizeDimension / cols;
    end

    % Resize the image
    resizedImage = imresize(image, scaleFactor);
end

% 
% 
% function image = readImage(imagePath)
%     % turn color/blackwhite image into [0, 255] double 2d array.
% 
%     % Construct the full file path
%     fullImagePath = fullfile(pwd, imagePath);
%     % fullImagePath
%     % Check if the file exists
%     if ~exist(fullImagePath, 'file')
%         error('File does not exist: %s', fullImagePath);
%     end
% 
%     % If the file exists, read it
%     rawImage = imread(fullImagePath);
%     image = double(rgb2grey1(rawImage));
% end
% 
% function gray = rgb2grey1(rgb)
%     if ndims(rgb) == 3
%         gray = 0.2989 * rgb(:,:,1) + 0.5870 * rgb(:,:,2) + 0.1140 * rgb(:,:,3);
%     else
%         gray = rgb;
%     end
% end
% 
    