% Define the range of k values and mu values
k_values = 1:10;
mu_a_values = [0.1, 0.25, 0.5];
mu_e_values = [0.5, 0.75, 1.0];
lambda = 0.5;
max_iter = 10000;
convergence_threshold = 1e-6;
a = 0.5; b = 2; c = 1;
priority_function = @(p, ta) p - a * (b * (1 - exp(-c * ta)));
% Initialize a matrix to store the W values for different k and mu values
W_values = zeros(length(k_values), length(mu_a_values));

% Loop over the k values and mu values
for k_idx = 1:length(k_values)
    for mu_idx = 1:length(mu_a_values)
        % Update the k, mu_a, and mu_e parameters
        k = k_values(k_idx);
        mu_a = mu_a_values(mu_idx);
        mu_e = mu_e_values(mu_idx);
        
        % Calculate the steady state probabilities P
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

        % Calculate Lq
        Lq = 0;
        for row = 1:k
            for i = 1:n
                Lq = Lq + (i + row) * P(i, row);
            end
        end

        % Calculate W
        Wq = Lq / (lambda);
        W = Wq + (1/mu_e);

        % Store the calculated W value in the matrix
        W_values(k_idx, mu_idx) = W;
    end
end

% Plot W vs k for each combination
figure;
hold on;
for c = 1:length(mu_a_values)
    plot(k_values, W_values(:, c), 'DisplayName', sprintf('mu_a=%d, mu_e=%d', mu_a_values(c), mu_e_values(c)));
end
hold off;

xlabel('k');
ylabel('W');
legend('show');
title('W vs k for different values of mu_a and mu_e with the same index');
