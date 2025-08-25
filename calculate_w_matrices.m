function [matrices, inverses] = calculate_w_matrices(NN, vD, am, w_max)
    % 输入:
    %   NN: 新增感染者数量向量 [N(t0), N(t0+1), ..., N(t0+max_index)]
    %   vD: 常数 ν + D
    %   a_plus_floor: 常数 ⌊a₊⌋ (通常为14)
    %   w_max: w的最大值
    % 输出:
    %   matrices: 元胞数组，存储每个w对应的矩阵
    %   inverses: 元胞数组，存储每个矩阵的逆矩阵

    % 矩阵大小为 (a_plus_floor + 1) × (a_plus_floor + 1)
    n = am + 1;
    matrices = cell(w_max, 1);
    inverses = cell(w_max, 1);
    
    for w = 1:w_max
        % 初始化矩阵
        A = zeros(n, n);
        
        % 构建前n-1列 (i = 1 到 n-1)
        for j=2:1:n
            threshold_col = n + 1 - j;
            for i=1:1:n
                if i<threshold_col
                    A(i,j)=NN(am+w-(j-1)+1);
                end
                if i==threshold_col
                    A(i,j)=2*NN(am+w-(j-1)+1);
                end
                if i>threshold_col && i<n
                    th1=i-threshold_col;
                    exp_val1 = th1 * vD;
                    I_val = ((1 - vD)^th1) * NN(am + w-(j-1)+1);
                    A(i,j) = I_val * exp(exp_val1);
                        
                    
                end

            end

        end
            
        
        for jj = 1:n
           
            A(n,jj)=((1-vD)^(jj-1))*(NN(am+w-(jj-1)+1))*(exp((jj-1)*vD));
            
        end
        
        % 存储矩阵
        %matrices{w} = A(2:end,2:end);%做个尝试
        matrices{w} = A;

        
        % 计算逆矩阵（添加条件数检查）
        %c = cond(A);
        %if c < 1e5
        %    inverses{w} = inv(A);
        %else
        %    warning('矩阵接近奇异 (w=%d, 条件数=%.2e)，使用伪逆', w, c);
        %    inverses{w} = pinv(A);
        %end
   % inverses{w} = pinv(A(2:end,2:end));%做个尝试
    inverses{w} = pinv(A);%做个尝试
    end
end