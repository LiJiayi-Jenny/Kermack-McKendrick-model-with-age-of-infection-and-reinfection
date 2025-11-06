function [matrices, inverses] =  calculate_inverse2(NN, am, I0)
    % inputs:
    %   NN: the vector of newly infected individuals [N(t0), N(t0+1), ..., N(t0+max_index)]
    %   vD: equals to constant v + D, v is the recovery rate of infecte
    %   indivudals, D is the mortality rate caused by the diseases
    %   a_plus_floor: costant ⌊a₊⌋ 
    %   w_max: the value of w
    % outputs:
    %   matrices: Cell array stores matrices of CoRm
    %   inverses: Cell array stores the inverses of matrices

    matrices = cell(am-1, 1);
    inverses = cell(am-1, 1);
    
    for m_val = 2:am
        % initialize matrix
        n = m_val+1;
        A = zeros(n, n);
        
        % construct the matrix
        for j=2:1:n
            threshold_col = n + 1 - j;
            for i=1:1:n-1
                if i<threshold_col
                    A(i,j)=NN(m_val-(j-1)+1);
                end
                if i==threshold_col
                    A(i,j)=2*NN(m_val-(j-1)+1);
                end
                if i>threshold_col && i<n
                    A(i,j)=NN(m_val-(j-1)+1);  
                end
            end

        end
        for jj = 1:n
           
            A(n,jj)=NN(m_val-(jj-1)+1);
            
        end
        
        % store the matrix
        matrices{m_val-1} = A;
        % caculate the inverses of matrix
        %c = cond(A);
        %if c < 1e5
        %    inverses{w} = inv(A);
        %else
        %    warning('The matrix is nearly singular (w=%d, condition number=%.2e), and the pseudoinverse is used', w, c);
        %    inverses{w} = pinv(A);
        %end
    inverses{m_val-1} = pinv(A);
    end
    newMatrix = [0, I0+NN(1);NN(2), I0];
        matrices = [{newMatrix}; matrices];
        innewMatrix=inv(newMatrix) ;
        inverses = [{innewMatrix}; inverses];
end
