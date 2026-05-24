clear;
clc;
close all;

% ------------------------------------------------------------
% Read convergence data
% ------------------------------------------------------------
data = readtable('CPP_Course\Exercises\Exercise2\output\convergence.csv');

N         = data.N;
h         = data.h;
max_error = data.max_error;
order     = data.order;
time_ms   = data.cpu_ms;

% Remove invalid rows, if any
valid = isfinite(h) & isfinite(max_error) & max_error > 0;

N_valid     = N(valid);
h_valid     = h(valid);
error_valid = max_error(valid);

% ------------------------------------------------------------
% Estimate slope from log-log data
% log(E) = p log(h) + c
% ------------------------------------------------------------
p = polyfit(log(h_valid), log(error_valid), 1);
estimated_order = p(1);

fprintf("Estimated convergence order from log-log fit: %.6f\n", estimated_order);

% ------------------------------------------------------------
% Reference O(h^2) line
% Scaled so it passes through the last error point
% ------------------------------------------------------------
C = error_valid(end) / h_valid(end)^2;
ref_error = C * h_valid.^2;

% ------------------------------------------------------------
% Figure 1: log-log convergence plot
% ------------------------------------------------------------
figure;

loglog(h_valid, error_valid, "o-", "LineWidth", 1.5, "MarkerSize", 7);
hold on;
loglog(h_valid, ref_error, "--", "LineWidth", 1.5);

grid on;
xlabel("h", "Interpreter", "latex");
ylabel("$E_\infty$", "Interpreter", "latex");
title("Convergence of the centered finite-difference method", ...
      "Interpreter", "latex");

legend( ...
    sprintf("Numerical error, slope = %.3f", estimated_order), ...
    "$O(h^2)$ reference", ...
    "Interpreter", "latex", ...
    "Location", "best" ...
);

set(gca, "FontSize", 12);

% Optional: reverse x-axis so refinement goes left to right visually
set(gca, "XDir", "reverse");

% Save figure
exportgraphics(gcf, "CPP_Course\Exercises\Exercise2/output/convergence_loglog.png", "Resolution", 300);

% ------------------------------------------------------------
% Figure 2: observed order
% ------------------------------------------------------------
valid_order = isfinite(order);

figure;

plot(N(valid_order), order(valid_order), "o-", ...
     "LineWidth", 1.5, "MarkerSize", 7);
hold on;
yline(2, "--", "Expected order 2", "LineWidth", 1.5);

grid on;
xlabel("$N$", "Interpreter", "latex");
ylabel("Observed order", "Interpreter", "latex");
title("Observed convergence order", "Interpreter", "latex");

set(gca, "FontSize", 12);

exportgraphics(gcf, "CPP_Course\Exercises\Exercise2/output/convergence_order.png", "Resolution", 300);

% ------------------------------------------------------------
% Figure 3: solve time
% ------------------------------------------------------------
figure;

plot(N, time_ms, "o-", "LineWidth", 1.5, "MarkerSize", 7);

grid on;
xlabel("$N$", "Interpreter", "latex");
ylabel("Time [ms]", "Interpreter", "latex");
title("Block Thomas solver runtime", "Interpreter", "latex");

set(gca, "FontSize", 12);

exportgraphics(gcf, "CPP_Course\Exercises\Exercise2/output/runtime.png", "Resolution", 300);