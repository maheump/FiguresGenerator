% This script plots entropy levels given p(A|B) and p(B|A) for observers
% learning either (i) transition probabilities, (ii) item frequency or
% (iii) alternation frequency.
% 
% Maxime Maheu, 11/2016

%% Initialization
%  ==============

% Clear the place
clear; close('all');

% Define grid
dt = 0.01;
t = 0:dt:1;
nt = numel(t);

% Create matrices of probability
pAgB = repmat(t,nt,1);
pBgA = repmat(t,nt,1)';

% Compute p(A) given p(A|B) and p(B|A)
pA = pAgB ./ (pAgB + pBgA);

% Compute p(alternation) given p(A|B) and p(B|A)
pAlt = (2.*pAgB.*pBgA) ./ (pAgB + pBgA);

% Define the entropy function
H = @(p) -(p .* log2(p) + (1-p) .* log2(1-p));

%% Display computed probabilities
%  ==============================

% Prepare the figure's window
figure('Position', [0 0.46 0.4 0.4]);

% Item frequency
subplot(1,2,1);
imagesc(t, t, pA); hold('on');
contour(t, t, pA, 'k', 'LineWidth', 1, 'ShowText', 'On');
title({'Item frequency',''});

% Alternation frequency
subplot(1,2,2);
imagesc(t, t, pAlt); hold('on');
contour(t, t, pAlt, 'k', 'LineWidth', 1, 'ShowText', 'On');
title({'Alternation frequency',''});

% Customize the...
for i = 1:2
    subplot(1,2,i);
    
    % ... axes
    cbr = colorbar('Location', 'SouthOutside'); cbr.LineWidth = 1;
    axis('square'); axis('xy'); caxis([0,1]);
    set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
    set(gca, 'FontSize', 15, 'LineWidth', 1);
    
    % ... labels
    cbr.Label.String = 'p';
    xlabel('p(A|B)'); ylabel('p(B|A)');
end
colormap(parula(1000));

%% Display the entropy
%  ===================

% Prepare the figure's window
figure('Position', [0.4 0.46 0.6 0.4]);

% Entropy in the space of transition probabilities
subplot(1,3,1);
imagesc(t, t, H(t)' * H(t)); hold('on');
contour(t, t, H(t)' * H(t), 'k', 'LineWidth', 1, 'ShowText', 'On');
title({'Transition probabilities',''});

% Entropy in the space of item frequency
subplot(1,3,2);
imagesc(t, t, H(pA)); hold('on');
contour(t, t, H(pA), 'k', 'LineWidth', 1, 'ShowText', 'On');
title({'Item frequency',''});

% Entropy in the space of alternation frequency
subplot(1,3,3);
imagesc(t, t, H(pAlt)); hold('on');
contour(t, t, H(pAlt), 'k', 'LineWidth', 1, 'ShowText', 'On');
title({'Alternation frequency',''});

% Customize the...
for i = 1:3
    subplot(1,3,i);
    
    % ... axes
    cbr = colorbar('Location', 'SouthOutside'); cbr.LineWidth = 1;
    axis('square'); axis('xy'); caxis([0,1]);
    set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
    set(gca, 'FontSize', 15, 'LineWidth', 1);
    
    % ... labels
    cbr.Label.String = 'H(p)';
    xlabel('p(A|B)'); ylabel('p(B|A)');
end
colormap(parula(1000));

%% Entropy on the diagonal
%  =======================

% Get the list of coefficients
w = 0:0.1:1;
nw = numel(w);

% Prepare the figure's window
figure('Position', [0 0.2 1 0.2]);

% For each calue of the coefficient
for i = 1:nw
    
    % Compute a weighted averaged between the two diagonals
    wp = (w(i).*pA) + (1-w(i)).*pAlt;
    
    % Plot the entropy of the resulting probability matrix
    subplot(1,nw,i);
    imagesc(t, t, H(wp));
    
    % Customize axes
    axis('square'); axis('xy'); caxis([0,1]);
    set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
    set(gca, 'FontSize', 15, 'LineWidth', 1);
    
    % Add labels to the plot
    if i == 1, xlabel('p(A|B)'); ylabel('p(B|A)'); end
    title(sprintf('w_{p(A)} = %1.2f', w(i)));
end
colormap(parula(1000));
