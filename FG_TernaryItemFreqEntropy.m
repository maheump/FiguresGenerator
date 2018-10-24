% This script displays Shannon entropy as a function of 3 orthogonal
% probability values
% 
% Maxime Maheu, 10/2018

% Define a probability grid
prec = 0.001;
grid = prec:prec:(1-prec);
[pA,pB] = meshgrid(grid, grid);

% Probabilities of A, B, and C must sum to 1
pC = 1 - (pA + pB);
idx = pC < 0;
pA(idx) = NaN;
pB(idx) = NaN;
pC(idx) = NaN;

% Concatenate probability values
p = cat(3, pA, pB, pC);
p = p ./ sum(p, 3, 'OmitNaN');

% Compute Shannon entropy on each set of probabilities
H = -sum(p .* log2(p), 3);

% Display the Shannon map
figure('Units', 'Normalized', 'Position', [0.4151 0.3058 0.1698 0.3500]);
surf(grid, grid, H, 'AlphaData', ~idx, 'EdgeColor', 'None', 'FaceColor', 'interp'); hold('on');
plot([1,0], [0,1], 'k-');
view([-45, 90]);
set(gca, 'XTick', [], 'YTick', []);
set(gca, 'Box', 'Off', 'Xgrid', 'Off', 'Ygrid', 'Off');
cbr = colorbar('Location', 'SouthOutside');
cbr.Label.String = 'Shannon entropy (bits)'; caxis([0, max(H(:))]);
xlabel('p(A)', 'FontSize', 15); ylabel('p(B)', 'FontSize', 15); ...
text(mean(grid), mean(grid), 'p(C)', 'FontSize', 15, 'VerticalAlignment', 'Bottom');