function [opt_path, conv_path] = Rosenbrock_gd_bisec(x0, y0, alpha_UB, epsilon, Nmax)
    xn = x0;
    yn = y0;
    opt_path = zeros(2, Nmax+1);
    conv_path = zeros(1, Nmax+1);
    for n = 1 : Nmax
        
        % call rosenbrock.m at (xn,yn) to get the current fn and gn
        [fn, gn] = rosenbrock(xn, yn);
        
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
        
        % normalize the gradient to get \hat{g}, the unit vector in the
        % gradient direction
        g_hat = gn / norm(gn);
        
        %% use bisection algorithm to find step size
        % define h'(alpha) for current iteration as described in main file
        hn_x =@(alpha) xn - alpha * g_hat(1);
        hn_y =@(alpha) yn - alpha * g_hat(2);
        
        g =@(alpha) [-2 * (1-hn_x(alpha)) - 4 * hn_x(alpha)*(hn_y(alpha)-hn_x(alpha).^2); 2 * (hn_y(alpha)-hn_x(alpha).^2)];
        
        h_prime =@(alpha) g(alpha)' * (-g_hat);
        
        % call bisection routine on h'(alpha) to find optimal step size
        [alphas, ~] = bisection(h_prime, 0, alpha_UB, 1e-10, 40);
        
        % take gradient descent step with optimal step size alpha. Note
        % that you must use the normalized gradient here because that's
        % what we used to calculate the optimal step size
        xn = xn - alphas(end) * g_hat(1); % replace with update
        yn = yn - alphas(end) * g_hat(2); % replace with update

    end
    opt_path(:, n+1) = [xn; yn];
    conv_path(:, n+1) = rosenbrock(xn, yn);
    opt_path = opt_path(:, 1:n+1);
    conv_path = conv_path(:, 1:n+1);
end
