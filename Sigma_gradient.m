function sigma_gradient = Sigma_gradient(z)
epsilon = 0.1;
sigma_gradient = z / sqrt(1 + epsilon * Sigma_norm(z));
end