% Define the range of k values and mu values
k_values = 1:10;
mu_a_values = [0.1, 0.25, 0.5];
mu_e_values = [0.5, 0.75, 1.0];
max_iter = 10000;
a = 0.5; b = 2; c = 1;
convergence_threshold = 1e-6;

priority_function = @(p, ta) p - a * (b * (1 - exp(-c * ta)));

% Initialize a matrix to store the gamma values for different k and mu values
gamma_values = zeros(length(k_values), length(mu_a_values));

% Loop over the k values and mu values
for k_idx = 1:length(k_values)
    for mu_idx = 1:length(mu_a_values)
        % Update the k, mu_a, and mu_e parameters
        k = k_values(k_idx);
        mu_a = mu_a_values(mu_idx);
        mu_e = mu_e_values(mu_idx);
        
        % Calculate the steady state probabilities P
        n = 100;
        lambda = 1/n;
        P = ones(n, k) / (n * k); % Initialize with uniform probabilities

        for iter = 1:max_iter
            P_old = P;

            for row = 1:k
                ta = row - 1;
                p_prime = priority_function(1:n, ta);

                for i = 1:n
                    sum_term = 0;
                    for j = 1:i
                        sum_term = sum_term + (i - j + 1) * lambda * p_prime(j) * P_old(j, row);
                    end

                    if i == 1
                        P(i, row) = (mu_a * p_prime(n) + mu_e * p_prime(n)) / sum_term;
                    else
                        P(i, row) = ((i - 1) * lambda * p_prime(i - 1)) / sum_term;
                    end
                end
            end

            P = P / sum(P(:));

            if max(abs(P(:) - P_old(:))) < convergence_threshold
                break;
            end
        end

        gamma = convergence_threshold/(mu_e/mu_a * P(n, k));

        % Store the calculated gamma value in the matrix
        gamma_values(k_idx, mu_idx) = gamma;
    end
end

% Plot gamma vs k for each combination
figure;
hold on;
for c = 1:length(mu_a_values)
    plot(k_values, gamma_values(:, c), 'DisplayName', sprintf('mu_a=%g', mu_a_values(c)));
end
hold off;

xlabel('k');
ylabel('gamma');
legend('show');
title('gamma vs k for different values of mu_a and mu_e');
