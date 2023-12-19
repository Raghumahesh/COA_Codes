% Parameters
lambda = 0.01;
k_values = 5:30;
r = 0.016;
mu_a_values = [0.1, 0.25, 0.5];
mu_e_values = [0.5, 0.75, 1.0];

n=5;
figure;
hold on;
grid on;

for idx = 1:length(mu_a_values)
    mu_a = mu_a_values(idx);
    mu_e = mu_e_values(idx);


    gamma_values = zeros(size(k_values));

    for k_idx = 1:length(k_values)
        k = k_values(k_idx);
        P = find_steady_state_probabilities(n, k, lambda, mu_a, mu_e);

        gamma = (mu_e/mu_a) * r;
        gamma_values(k_idx) = gamma;
    end

    % Plot the graph for the current mu_a and mu_e values
    plot(k_values, gamma_values, '-o', 'DisplayName', ['\mu_a = ' num2str(mu_a) ', \mu_e = ' num2str(mu_e)]);
end

xlabel('k');
ylabel('Gamma');
title('Gamma vs. k for different values of \mu_a and \mu_e');
legend('show');
hold off;


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