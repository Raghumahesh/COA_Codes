% Parameters
lambda = 0.1;
n = 10;

% Define the range of k values and mu values
k_values = 1:10;
mu_a_values = [0.1, 0.25, 0.5];
mu_e_values = [0.5, 0.75, 1.0];
% Initialize a matrix to store the Wq values for different k and mu values
Wq_values = zeros(length(k_values), length(mu_a_values));

% Loop over the k values and mu values
for k_idx = 1:length(k_values)
    for mu_idx = 1:length(mu_a_values)
        % Update the k, mu_a, and mu_e parameters
        k = k_values(k_idx);
        mu_a = mu_a_values(mu_idx);
        mu_e = mu_e_values(mu_idx);

        % Calculate the steady state probabilities P
        P = find_steady_state_probabilities(n, k, lambda, mu_a, mu_e);

        % Calculate Lq
        Lq1 = 0;
        for row = 1:n+1
            for col = 1:k+1
                Lq = Lq + (row - 1 + col - 1) * P(row, col);
                Wq_values(k_idx, mu_idx) = Lq;
            end
        end
    end
end

% Plot Wq vs k for each combination
figure;
hold on;
for c = 1:length(mu_a_values)
    plot(k_values, Wq_values(:, c), 'DisplayName', sprintf('mu_a=%d, mu_e=%d', mu_a_values(c), mu_e_values(c)));
end
hold off;

xlabel('k');
ylabel('Lq');
legend('show');
title('Lq vs k for different values of mu_a and mu_e with the same index');

function [P] = find_steady_state_probabilities(n, k, lambda, mu_a, mu_e)
    num_states = (n+1)*(k+1);
    A = zeros(num_states);

    % Populate matrix A
    for row = 0:n
        for col = 0:k
            state_idx = row*(k+1) + col + 1;
            if row == n && col == k
                A(state_idx, state_idx) = 1;
            else
                if row < n
                    A(state_idx, state_idx) = -(row+1)*lambda;
                    if row > 0
                        A(state_idx, state_idx - (k+1)) = row*lambda;
                    end
                end
                if col < k
                    A(state_idx, state_idx) = A(state_idx, state_idx) - mu_e;
                    A(state_idx, state_idx + 1) = mu_e;
                end
                if col > 0
                    A(state_idx, state_idx) = A(state_idx, state_idx) - mu_a;
                    A(state_idx, state_idx - 1) = mu_a;
                end
            end
        end
    end

    % Add normalization equation
    A(end, :) = ones(1, num_states);
    b = zeros(num_states, 1);
    b(end) = 1;

    % Solve for P
    P = A \ b;
    P = reshape(P, [n+1, k+1]);
end
