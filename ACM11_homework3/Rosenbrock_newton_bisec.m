function [opt_path, conv_path] = Rosenbrock_newton_bisec(x0, y0, alpha_UB, epsilon, Nmax)
    xn = x0;
    yn = y0;
    opt_path = zeros(2, Nmax+1);
    conv_path = zeros(1, Nmax+1);
    for n = 1 : Nmax
        % call rosenbrock.m at (xn,yn) to get the current fn, gn, and Hn
        [fn, gn, Hn] = rosenbrock(xn, yn);
        
        % store fn in nth entry of conv_path, and store xn,yn in nth column
        % of opt path
        conv_path(1, n) = fn;
        opt_path(1, n) = xn;
        opt_path(2, n) = yn;
        
        % check to see if norm of gradient is below epsilon tolerance - if
        % so, break the for loop.
        if norm(gn) < epsilon
            break;
        end        
        
        %% use bisection algorithm to find step size
        % define \hat{g}, the Newton update direction, as defined in the
        % main file - note that this is NOT the same as \hat{g} for
        % gradient descent
        dn = Hn \ gn;
        g_hat = dn / norm(dn);

        % check if \hat{g} is a descent direction. If not, change the
        % sign of \hat{g}
        if (gn') * g_hat < 0
            g_hat = (-1) * g_hat;
        end 
        
        % define q(alpha) for current iteration as described in main file
        q_x =@(alpha) xn - alpha * g_hat(1);
        q_y =@(alpha) yn - alpha * g_hat(2);
        g =@(alpha) [-2 * (1-q_x(alpha)) - 4 * q_x(alpha)*(q_y(alpha)-q_x(alpha).^2); 2 * (q_y(alpha)-q_x(alpha).^2)];

        q =@(alpha) g(alpha)' * (-1)*g_hat;
        
        % call bisection routine on q(alpha) to find optimal step size
        [alphas, ~] = bisection(q, 0, alpha_UB, 1e-10, 40);
        
        % take Newton step with optimal step size alpha (you must use the
        % normalized Newton direction here)
        xn = xn - alphas(end) * g_hat(1); % replace with update
        yn = yn - alphas(end) * g_hat(2); % replace with update

    end
    opt_path(:, n+1) = [xn; yn];
    conv_path(:, n+1) = rosenbrock(xn,yn);
    opt_path = opt_path(:, 1:n+1);
    conv_path = conv_path(:, 1:n+1);
end