function W = poly_with_chebyshev(z, max_order)
    Sx = 512;
    Sy = 512;
    
    % 归一化 z 的实部和虚部到 [-1, 1] 区间
    x = real(z);
    y = imag(z);
    x_norm = 2 * (x - min(x(:))) / (max(x(:)) - min(x(:))) - 1;  % 映射到 [-1, 1]
    y_norm = 2 * (y - min(y(:))) / (max(y(:)) - min(y(:))) - 1;  % 映射到 [-1, 1]

    % 初始化 W 为合适的大小
    % 这里的basis函数的数量是 (max_order+1)*(max_order+2)/2
    W = zeros([size(z, 1), (max_order + 1) * (max_order + 2) / 2]);
    id = 1;
    
    % 计算基函数顺序为 (m, n) 的对称排列
    for max_t = 0:max_order
        for m = 0:max_t
            n = max_t - m;  % 保证 m + n = max_t

            % 计算第 m 和 n 阶的切比雪夫多项式 T_m(x) 和 T_n(y)
            Tm_x = chebyshev_polynomial(x_norm, m);
            Tn_y = chebyshev_polynomial(y_norm, n);
            
            % 计算基函数的乘积
            W(:, id) = Tm_x .* Tn_y;

            % 保持基函数不归一化，以保持正交性
            % 归一化到 [0, Sx] 和 [0, Sy] 范围仅在最后
            if id > 1
                W(:, id) = (W(:, id) - min(W(:, id))) / (max(W(:, id)) - min(W(:, id))) * Sx;  % 对x轴进行归一化
                W(:, id) = (W(:, id) - min(W(:, id))) / (max(W(:, id)) - min(W(:, id))) * Sy;  % 对y轴进行归一化
            end

            id = id + 1;
        end
    end
end

% 计算切比雪夫多项式
function T = chebyshev_polynomial(x, n)
    if n == 0
        T = ones(size(x));  % T_0(x) = 1
    elseif n == 1
        T = x;  % T_1(x) = x
    else
        T_0 = ones(size(x));
        T_1 = x;
        for k = 2:n
            T = 2 * x .* T_1 - T_0;  % T_{n+1}(x) = 2x T_n(x) - T_{n-1}(x)
            T_0 = T_1;
            T_1 = T;
        end
    end
end
