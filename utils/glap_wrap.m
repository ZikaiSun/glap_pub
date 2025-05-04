
function [I1_wrap, I2, mask, u_total] = glap_wrap(I1_original, I2_original, sigma_order_list)


% Initialize separate lists for sigma and order
sigma_list = [];
order_list = [];
for i = 1:length(sigma_order_list)
    pair = sigma_order_list{i};
    sigma_list = [sigma_list, pair(1)];
    order_list = [order_list, pair(2)];
end
max_order = max(order_list);
iteration  = length(sigma_list);

    resize_list = ones(size(sigma_list))*512;
    scale_base  = 512; 
    sigma_base = 32; 
    k = 0.5; 
    num_sigma = 1;
    sigma_mask_list = [sigma_base*k^3, sigma_base*k^2, sigma_base*k, sigma_base];
    sigma_mask_list = sigma_mask_list((end-num_sigma):end);
    [M_original, N_original] = size(I1_original);
    mask = ones([M_original, N_original]);
    coeffs = zeros([(max_order +1)*(max_order +2)/2,1]);
    scale_i = resize_list(1); 
    I1 = I1_original;
    I2 = I2_original;
    mask = imresize(mask,[scale_i scale_i]);
    [M, N]=size(I1);
    z = repmat((1:M)',1,N)+1i*repmat((1:N),M,1);

    for i=1:iteration
        mask_common = shrink(logical(imresize(ones([M, N]),[M, N])), round(sigma_list(i)*2* (sqrt(M*N)) /scale_base));
        W = poly_with_chebyshev(z(:), order_list(i));
        [~, wn] = size(W);
        u_total = reshape(-4*W*coeffs(1:wn), [M,N]);
        I1_wrap = imshift(I1, -u_total,'bilinear'); 
        maskrw = imshift(mask_common, -u_total,'nearestneighbor').*mask_common > 0.01;
        I1_intensity_wrap = histomatch(I1_wrap, I2, imshift(ones(size(I2)), -u_total,'nearestneighbor') > 0.01);

        % 2: occlusion estimation
        mask_sum = zeros(size(I1_wrap));
        for sigma_i = sigma_mask_list
            [A, b, c] = get_Abc(I1_intensity_wrap, I2, maskrw, sigma_i); % this function only use the first 3 parameters coefficients.
            err = reshape(abs(A(:,:)*c-b(:,:)), size(I1_wrap));
            mask_i = err <  3* median(err(:));  % err <  0.3 * max(err(:));
            mask_sum = mask_sum + mask_i;
        end
        mask = ((mask_sum >= (size(sigma_mask_list,2)-1)) .*  maskrw) > 0.001;

        % 3: GLAP algorithm
        [~, c] = glap_a(I1_intensity_wrap, I2, mask, sigma_list(i)*scale_i/scale_base, W, max_order);
        coeffs = coeffs + c;

    end
    u_total = reshape(-4*W*coeffs(1:wn), [M,N]);
    I1_wrap = imshift(I1, -u_total,'cubicOMOMS'); 
end


function [u, c] = glap_a(I1, I2, mask, sigma, W, order)
    mask_index = logical(mask(:));

    Gauss_LP=@(wx,wy)exp(-sigma^2*(wx.^2+wy.^2)/2);
    Gauss_HP=@(wx,wy)1i*(wx-1i*wy).*exp(-sigma^2*(wx.^2+wy.^2)/2);
    
    sublow = imagefilter(I2-I1,Gauss_LP); 
    addhigh = imagefilter(I2+I1,Gauss_HP);
    b=-reshape(sublow,[],1);
    A0=reshape(addhigh,[],1);
    A1=A0.*W;
    A1=[A1,conj(A1)];
    warning('off','all')
    coeffs = A1(mask,:)\b(mask);
    [~, wn] = size(W);
    u = 0;
    c = zeros([(order+1)*(order+2)/2,1]);
    c(1:wn) = coeffs(1:wn);
end
        
function [A, b, coeffs] = get_Abc(I1, I2, mask, sigma)
    [M,N]=size(I1);
    Gauss_LP=@(wx,wy)exp(-sigma^2*(wx.^2+wy.^2)/2);
    Gauss_HP=@(wx,wy)1i*(wx-1i*wy).*exp(-sigma^2*(wx.^2+wy.^2)/2);
    sublow = imagefilter(I2-I1,Gauss_LP); 
    addhigh = imagefilter(I2+I1,Gauss_HP);
    b=-reshape(sublow,[],1);
    A0=reshape(addhigh,[],1);
    z=repmat((1:M)',1,N)+1i*repmat((1:N),M,1);
    A=A0.*[ones(M*N,1),z(:),conj(z(:))];
    A=[A,conj(A)];
    coeffs = A(mask,:)\b(mask);
end
    
    