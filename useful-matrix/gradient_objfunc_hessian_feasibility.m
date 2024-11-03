% Compute gradient and hessian
function [grad, hess] = compute_grad_hess(H, g, A, b, u, t)
    d = 1./(b - A*u);
    grad = t*(H*u + g) + A'*d;
    hess = t*H + A' * diag(d.^2) * A;
end

% Objective function
function J = compute_J(H, g, A, b, u, t)
    J = t*(0.5*u'*H*u + g'*u) - sum(log(b - A*u));
end

% Verify if is in a feasible domain
function in_domain = is_in_domain(A, b, u, descent, t_newton)
    in_domain = all(A*(u + t_newton*descent) < b);
end
