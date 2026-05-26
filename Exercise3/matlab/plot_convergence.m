clear;
clc;
close all;

% ============================================================
% Comparison of Copy, View, and Map implementations
% Block tridiagonal Thomas-LDU solver
% ============================================================

base_path = "CPP_Course\Exercises\Exercise3\output\";

files = [
    "convergence_copy.csv"
    "convergence_view.csv"
    "convergence_map.csv"
];

labels = [
    "Copy"
    "View"
    "Map"
];

% ============================================================
% Figure 1: log-log convergence error
% Maximum error vs grid size h
% ============================================================

figure;
hold on;

h_ref = [];
error_ref = [];

for k = 1:numel(files)

    data = readtable(base_path + files(k));

    h         = data.h;
    max_error = data.max_error;

    valid = isfinite(h) & isfinite(max_error) & max_error > 0;

    h_valid     = h(valid);
    error_valid = max_error(valid);

    % Sort by h so the curve is connected correctly
    [h_valid, idx] = sort(h_valid, "descend");
    error_valid = error_valid(idx);

    % Estimate slope from log-log data
    p = polyfit(log(h_valid), log(error_valid), 1);
    estimated_order = p(1);

    fprintf("%s estimated convergence order: %.6f\n", ...
            labels(k), estimated_order);

    plot(h_valid, error_valid, "o-", ...
         "LineWidth", 1.5, ...
         "MarkerSize", 7, ...
         "DisplayName", sprintf("%s, slope = %.3f", labels(k), estimated_order));

    if isempty(h_ref)
        h_ref = h_valid;
        error_ref = error_valid;
    end
end

% Reference O(h^2) line
C = error_ref(end) / h_ref(end)^2;
ref_error = C * h_ref.^2;

plot(h_ref, ref_error, "--", ...
     "LineWidth", 1.5, ...
     "DisplayName", "$O(h^2)$ reference");

grid on;
xlabel("h", "Interpreter", "latex");
ylabel("$E_\infty$", "Interpreter", "latex");
title("Convergence of the block tridiagonal Thomas-LDU solver", ...
      "Interpreter", "latex");

legend("Interpreter", "latex", "Location", "best");

set(gca, "FontSize", 12);

% Force true log-log axes
set(gca, "XScale", "log");
set(gca, "YScale", "log");

% Refinement goes left to right visually
set(gca, "XDir", "reverse");

exportgraphics(gcf, ...
    base_path + "comparison_convergence_loglog.png", ...
    "Resolution", 300);

% ============================================================
% Figure 2: observed convergence order
% ============================================================

figure;
hold on;

for k = 1:numel(files)
    data = readtable(base_path + files(k));

    N     = data.N;
    order = data.order;

    valid_order = isfinite(order);

    plot(N(valid_order), order(valid_order), "o-", ...
         "LineWidth", 1.5, ...
         "MarkerSize", 7, ...
         "DisplayName", labels(k));
end

yline(2, "--", "Expected order 2", ...
      "LineWidth", 1.5, DisplayName='ref');

grid on;
xlabel("$N$", "Interpreter", "latex");
ylabel("Observed order", "Interpreter", "latex");
title("Observed convergence order", "Interpreter", "latex");

legend("Interpreter", "latex", "Location", "best");
set(gca, "FontSize", 12);

exportgraphics(gcf, ...
    base_path + "comparison_convergence_order.png", ...
    "Resolution", 300);

% ============================================================
% Figure 3: CPU time
% ============================================================

figure;
hold on;

for k = 1:numel(files)

    data = readtable(base_path + files(k));

    N = data.N;
    time_ms = get_cpu_time(data);

    valid_time = isfinite(N) & isfinite(time_ms) & time_ms > 0;

    plot(N(valid_time), time_ms(valid_time), "o-", ...
         "LineWidth", 1.5, ...
         "MarkerSize", 7, ...
         "DisplayName", labels(k));
end

grid on;
xlabel("$N$", "Interpreter", "latex");
ylabel("CPU time [ms]", "Interpreter", "latex");
title("Block Thomas-LDU solver runtime", "Interpreter", "latex");

legend("Interpreter", "latex", "Location", "best");
set(gca, "FontSize", 12);

exportgraphics(gcf, ...
    base_path + "comparison_runtime.png", ...
    "Resolution", 300);

% ============================================================
% Local function
% ============================================================

function time_ms = get_cpu_time(data)

    names = data.Properties.VariableNames;

    if ismember("cpu_ms_avg", names)
        time_ms = data.cpu_ms_avg;
    elseif ismember("cpu_ms", names)
        time_ms = data.cpu_ms;
    elseif ismember("time_ms_avg", names)
        time_ms = data.time_ms_avg;
    elseif ismember("time_ms", names)
        time_ms = data.time_ms;
    else
        error("No CPU time column found. Expected cpu_ms, cpu_ms_avg, time_ms, or time_ms_avg.");
    end

end