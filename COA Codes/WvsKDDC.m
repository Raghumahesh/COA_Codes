% Parameters
lambda = 0.5;


k_values = 1:10;
mu_a_values = [0.1, 0.25, 0.5];
mu_e_values = [0.5, 0.75, 1.0];

figure;
hold on;
grid on;

for mu_a_idx = 1:length(mu_a_values)
    mu_a = mu_a_values(mu_a_idx);
    mu_e = mu_e_values(mu_a_idx);

    W_values = zeros(size(k_values));

    for idx = 1:length(k_values)
        k = k_values(idx);
        P = find_steady_state_probabilities(n, k, lambda, mu_a, mu_e);

        L = 0;
        Lq = 0;
        W=0;
        Wq=0;
        for row = 0:n
            for col = 0:k
                Lq = Lq + (row + col) * P(row+1, col+1);
                Wq= Lq / lambda;
                W = W + (1/mu_e);
            end
        end

        W_values(idx) = W;
    end

    % Plot the graph for the current mu_a and mu_e values
    plot(k_values, W_values, '-o', 'DisplayName', ['\mu_a = ' num2str(mu_a) ', \mu_e = ' num2str(mu_e)]);
end

xlabel('k');
ylabel('W');
title('W vs. k for different values of \mu_a and \mu_e');
legend('show');
hold off;


% Calculate L, Lq, W, Wq, and gamma


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