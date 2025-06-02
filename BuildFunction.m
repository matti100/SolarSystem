function F = BuildFunction(x, M, G)

    % Initialize right-hand side function
    [n, ~] = size(x);
    N = n / 6; % number of bodies
    F = zeros(n, 1);
    
    for n = 1:N
        i = 6*(n-1) + 1;

        r_i = x(i:i+2);
        v_i = x(i+3:i+5);
    
        % r_i_dot
        F(i:i+2) = v_i;
    
        % v_i_dot
        f = zeros(3,1);
        for m = 1:N
            
            if m ~= n
                j = 6*(m-1) + 1;

                r_j = x(j:j+2);
                r_ji = r_i - r_j;

                d = norm(r_ji) + eps;

                f = f - M(m) .* r_ji ./ (d^3);
            end
    
        end
    
        F(i+3:i+5) = G .* f;
    
    end

end