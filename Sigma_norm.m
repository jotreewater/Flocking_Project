function sigma_norm = Sigma_norm(z)
epsilon = 0.1;
sigma_norm = (1 / epsilon) * (sqrt(1 + epsilon *(norm(z))^2) - 1);
end