clear;
clc;
close all;

files = dir("../output/*.csv");

Nmax = 300;

results = table();

for k = 1:length(files)

    filepath = fullfile(files(k).folder, files(k).name);
    T = readtable(filepath);

    N = min(Nmax, height(T));

    r = T.rank(1:N);
    f = T.frequency(1:N);

    % Normalize frequency
    f = f / sum(f);

    % Zipf-Mandelbrot model:
    % f(r) = c / (r + beta)^alpha
    model = @(p, r) p(1) ./ ((r + p(3)).^p(2));

    % Least squares error
    error_fun = @(p) sum((f - model(p, r)).^2);

    % Initial guess:
    % p(1) = c
    % p(2) = alpha
    % p(3) = beta
    p0 = [f(1), 1, 2.7];

    p = fminsearch(error_fun, p0);

    c     = p(1);
    alpha = p(2);
    beta  = p(3);

    f_fit = model(p, r);

    name = erase(files(k).name, ".csv");

    fprintf("\n%s\n", name);
    fprintf("c     = %.6e\n", c);
    fprintf("alpha = %.6f\n", alpha);
    fprintf("beta  = %.6f\n", beta);

    results = [results;
        table(string(name), c, alpha, beta, ...
        'VariableNames', {'File','c','alpha','beta'})];

    figure;

    loglog(r, f, "o");
    hold on;
    loglog(r, f_fit, "-", "LineWidth", 2);

    grid on;
    xlabel("Rank");
    ylabel("Normalized Frequency");
    title("Zipf-Mandelbrot Fit - " + string(name));
    legend("Data", "Fit", "Location", "best");

    txt = sprintf(['f(r) = c / (r + \\beta)^\\alpha\n' ...
                   'Least squares\n' ...
                   'c = %.3e\n' ...
                   '\\alpha = %.3f\n' ...
                   '\\beta = %.3f'], ...
                   c, alpha, beta);

    text(0.05, 0.05, txt, ...
         'Units', 'normalized', ...
         'FontSize', 11, ...
         'BackgroundColor', 'white');

    pdf_name = "../output/" + string(name) + "_fit.pdf";

    exportgraphics(gcf, pdf_name, ...
                   "ContentType", "vector");

    fprintf("Saved figure: %s\n", pdf_name);

    hold off;
end

disp(results);