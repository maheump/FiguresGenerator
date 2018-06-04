% Script illustrating the concept of model evidence used in Bayesian Model
% Selection
% 
% Maxime Maheu, 2016

% Precision of the grid
xgrid = 0:0.01:1;
subgrid = round(linspace(1,numel(xgrid),40));

% Widths
lw1 = 4;
lw2 = 2;

% Sizes
fontsize = 15;
marksize = 8;

% Colors
cols  = lines(2);
lgray = repmat(200 / 255, 1, 3);
dgray = repmat(128 / 255, 1, 3);

% Create distributions that will be used as likelihood proxi
y1 = betapdf(xgrid, 5, 5);
y2 = betapdf(xgrid, 5, 5);
y3 = repmat(y1', 1, numel(y2)) .* repmat(y2, numel(y1), 1);

% Normalize distributions
y3 = y3 ./ sum(y3(:));
y2 = y2 ./ (sum(y2(:)).*60);

figure; lgd = NaN(1,2);

% Plot the likelihood of the complex model
lgd(1) = surf(xgrid(subgrid), xgrid(subgrid), y3(subgrid,subgrid), 'EdgeColor', 'k');
hold('on'); colormap(gray);

% Plot the likelihood of the simple model
plot3(xgrid, zeros(1,numel(y2)), y2, 'r-', 'LineWidth', lw2);
lgd(2) = fill3(xgrid, zeros(1,numel(y2)), y2, 'k', ...
    'EdgeColor', 'None', 'FaceColor', 'r', 'FaceAlpha', 0.5);

% Plot the maxima and their values
plot3(0.5, 0.5, max(y3(:)), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 15);
text(0.5, 0.5, max(y3(:)), sprintf('  Mean = %1.5f\n   Max = %1.5f', mean(y3(:)), ...
    max(y3(:))), 'HorizontalAlignment', 'Left', 'FontSize', 15, 'Color', 'b');
text(0.5, 0, max(y2(:)), sprintf('  Mean = %1.5f\n   Max = %1.5f', mean(y2(:)), ...
    max(y2(:))), 'HorizontalAlignment', 'Left', 'FontSize', 15, 'Color', 'b');
plot3(0.5, 0, max(y2(:)), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 15);

% Customize the axes
axis([0,1,0,1,0,max(y3(:))]);
set(gca, 'XTickLabel', {}, 'YTickLabel', {}, 'ZTickLabel', {}, ...
    'Xgrid', 'Off', 'Ygrid', 'Off', 'FontSize', 15, 'LineWidth', 1);

% Add labels
legend(lgd, {'Model 1', 'Model 2'}, 'Location', 'NorthEast');
xlabel('\theta_{1}');
ylabel('\theta_{2}');
zlabel('p( y | M, \theta_{1}, \theta_{2} )');