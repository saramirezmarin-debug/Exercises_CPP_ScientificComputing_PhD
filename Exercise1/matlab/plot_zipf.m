% Get all csv files
files = dir('../output/*.csv');

% Create figure
figure;

set(gcf, 'Color', 'w');
set(gcf, 'Position', [100 100 900 700]);

hold on;
grid on;
box on;

% Axis properties
ax = gca;

ax.FontSize = 14;
ax.LineWidth = 1.5;

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

% Automatic colors
colors = lines(length(files));

for k = 1:length(files)

    filepath = fullfile(files(k).folder, ...
                        files(k).name);

    T = readtable(filepath);

    % Rank
    r = T.rank;

    % Normalized frequency
    f = T.frequency;
    f = f / sum(f);

    % Remove .csv extension
    name = erase(files(k).name, '.csv');

    % Plot Zipf curve
    loglog(r, ...
           f, ...
           'Color', colors(k,:), ...
           'LineWidth', 2, ...
           'DisplayName', name);

end

xlabel('Rank (log scale)', ...
       'FontSize', 16);

ylabel('Normalized Frequency (log scale)', ...
       'FontSize', 16);

title('Zipf', ...
      'FontSize', 18);

legend('Location', 'southwest', ...
       'Interpreter', 'none', ...
       'FontSize', 12);

xlim([1 inf]);

hold off;

% Export PDF
exportgraphics(gcf, ...
               '../output/all_zipf.pdf', ...
               'ContentType', 'vector');

fprintf('Saved figure: all_zipf.pdf\n');