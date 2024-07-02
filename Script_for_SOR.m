A = [
    -2.011, 1, 0, 0, 0, 0, 0, 0, 0;
    1, -2.012, 1, 0, 0, 0, 0, 0, 0;
    0, 1, -2.013, 1, 0, 0, 0, 0, 0;
    0, 0, 1, -2.014, 1, 0, 0, 0, 0;
    0, 0, 0, 1, -2.015, 1, 0, 0, 0;
    0, 0, 0, 0, 1, -2.016, 1, 0, 0;
    0, 0, 0, 0, 0, 1, -2.017, 1, 0;
    0, 0, 0, 0, 0, 0, 1, -2.018, 1;
    0, 0, 0, 0, 0, 0, 0, 1, -2.019
];

b = [
    -0.994974;
    1.57407e-3;
    -8.96677e-4;
    -2.71137e-3;
    -4.07407e-3;
    -5.11719e-3;
    -5.92917e-3;
    -6.57065e-3;
    -0.507084
];

x0 = [
    0.95;
    0.9;
    0.85;
    0.8;
    0.75;
    0.7;
    0.65;
    0.6;
    0.55
];

d = 1; % Band width
tol = 1e-4; 
nmax = 100; % Maximum # of iterations

w_values = 1.0:0.1:1.9; % Weighs
iterations = zeros(size(w_values));

for i = 1:length(w_values)
    w = w_values(i);
    [~, nit] = sor(A, b, x0, w, d, tol, nmax);
    iterations(i) = nit;
end

figure;
plot(w_values, iterations, '-*');
xlabel('Relaxation parameter w');
ylabel('Number of iterations');
title('Number of iterations Vs Relaxation w');
grid on;

[~, min_idx] = min(iterations);
optimal_w = w_values(min_idx); % Best weight
fprintf('Optimal relaxation parameter w is = %.2f\n', optimal_w);