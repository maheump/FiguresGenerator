% This script displays Shannon entropy as a function of 3 orthogonal
% probability values
% 
% Maxime Maheu, 10/2018

% Define a probability grid
prec = 0.01;
grid = prec:prec:(1-prec);
[pA,pB,pC] = meshgrid(grid, grid, grid);

% Compute Shannon entropy on each set of probabilities
p = cat(4, pA, pB, pC);
H = -sum(p .* log2(p), 4);

% Display the Shannon map
figure;
h = slice(grid, grid, grid, H, grid, grid, grid);
set(h, 'EdgeColor','None', 'FaceColor','interp');
alpha(0.05);
axis('equal');
xlabel('p(A)'); ylabel('p(B)'); zlabel('p(C)');
cbr = colorbar; cbr.Label.String = 'Shannon entropy (bits)'; caxis([0, max(H(:))]); 