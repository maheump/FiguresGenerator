% Script generating small cartoons illustrating which quantities one can
% derive from a Bayesian inference on a sequence.
% 
% Maxime Maheu, 11/2016

%% Initialization
%  ==============

% Precision of the grid
xgrid = 0:0.01:1;

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

%
figure('Name', 'Cartoons for quantities', 'Position', [0 0.0425 0.1005 0.8775]);

% Blank subplot
subplot(9,1,1);
axis('off');

%% Mean
%  ====

subplot(9,1,2);
y1 = betapdf(xgrid, 3, 7);
m = mean(y1.*xgrid);
plot(repmat(m,1,2), [0,max(y1)], 'k:', 'LineWidth', lw2); hold('on');
plot(xgrid, y1, 'Color', dgray, 'LineWidth', lw1);
ylim([0,max(y1)]);
ylabel('Mean');

%% Mode
%  ====

subplot(9,1,3);
y1 = betapdf(xgrid, 3, 7);
[m,idx1] = max(y1);
plot(xgrid, y1, 'Color', dgray, 'LineWidth', lw1); hold('on');
plot(xgrid(idx1), m, 'ko', 'LineWidth', lw2, 'MarkerFaceColor', 'w', 'MarkerSize', marksize);
axis([0,1,0,max(y1)]);
ylabel('Mode');

%% Variance
%  ========

subplot(9,1,4);
y1 = betapdf(xgrid, 5, 5);
y2 = betapdf(xgrid, 50, 50);
plot(xgrid, y1, 'LineWidth', lw1); hold('on');
plot(xgrid, y2, 'LineWidth', lw1); 
plot([mean(xgrid.*y1)-std(xgrid.*y1), mean(xgrid.*y1)+std(xgrid.*y1)], repmat(max(y2)+1,1,2), '-', 'LineWidth', lw1);
plot([mean(xgrid.*y2)-std(xgrid.*y2), mean(xgrid.*y2)+std(xgrid.*y2)], repmat(max(y2)+1,1,2), '-', 'LineWidth', lw1);
ylim([0,max(y2)+1]);
ylabel('Precision');

%% Prediction
%  ==========

subplot(9,1,5);
y1 = betapdf(xgrid, 2, 4);
y2 = betapdf(xgrid, 20, 4);
y3 = betapdf(xgrid, 2+20, 4+4);
m = nanmean(y3.*xgrid);
plot(xgrid, y1, ':', 'LineWidth', 2); hold('on');
plot(xgrid, y2, ':', 'LineWidth', 2); 
plot(repmat(m,1,2), [0,max(y3)], 'k:', 'LineWidth', lw2);
plot(xgrid, y3, '-', 'LineWidth', lw1); 
ylim([0,max(y2)]);
ylabel('Prediction');

%% Entropy
%  =======

subplot(9,1,6);
H = @(p) -(p .* log2(p) + (1-p) .* log2(1-p));
y1 = H(xgrid(2:end-1));
[~,idx1] = min(abs(xgrid-m));
plot(xgrid(2:end-1), y1, 'LineWidth', lw1); hold('on');
plot(xgrid(idx1), y1(idx1), 'ko', 'LineWidth', lw2, 'MarkerFaceColor', 'w', 'MarkerSize', marksize);
ylim([0,max(y1)]);
ylabel('Entropy');

%% Surprise
%  ========

subplot(9,1,7);
event = 60;
plot([0,xgrid(event)], repmat(y3(event),1,2), 'k:', 'LineWidth', lw2); hold('on');
plot(xgrid, y3, 'LineWidth', lw1); hold('on');
plot(repmat(xgrid(event),1,2), [max(y3), max(y3)/2], 'k-', 'LineWidth', lw2);
plot(xgrid(event), max(y3)/2, 'kv', 'MarkerFaceColor', 'k');
plot(xgrid(event), y3(event), 'ko', 'LineWidth', lw2, 'MarkerFaceColor', 'w', 'MarkerSize', marksize);
ylabel('Surprise');

%% Update
%  ======

subplot(9,1,8);
y1 = betapdf(xgrid, 3, 2);
y2 = betapdf(xgrid, 3, 3);
[~,idx1] = max(y1);
[~,idx2] = max(y2);
[~,idx3] = min(abs(y2-y1));
idx3 = 61;
lred = [255 150 150] ./ 255;
dred = [200 000 000] ./ 255;
fill([xgrid, fliplr(xgrid)], [y1(1:idx3),y2(idx3+1:end),y1(end:-1:idx3),y2(idx3-1:-1:1)], ...
    'k', 'EdgeColor', 'None', 'FaceColor', mean([lred;dred],1), 'FaceAlpha', 0.3); hold('on');
plot(xgrid, y1, 'LineWidth', lw1, 'Color', lred);
plot(xgrid, y2, 'LineWidth', lw1, 'Color', dred);
plot(xgrid([idx1,idx2]), repmat(max(y2)+0.5,1,2), 'k-', 'LineWidth', lw2);
plot(xgrid(idx2), max(y2)+0.5, 'k<', 'MarkerFaceColor', 'k');
%plot(grid, y2-y1, 'LineWidth', lw1); 
ylim([0,max(y2)+0.5]);
ylabel('Update');

%% Model evidence
%  ==============

subplot(9,1,9);
subgrid = round(linspace(1,numel(xgrid),21));
y1 = betapdf(xgrid, 5, 5);
y2 = betapdf(xgrid, 5, 5);
y3 = repmat(y1', 1, numel(y2)) .* repmat(y2, numel(y1), 1);
surf(xgrid(subgrid), xgrid(subgrid), y3(subgrid,subgrid), 'EdgeColor', 'k'); hold('on'); colormap(flipud(lgray));
plot3(xgrid, zeros(1,numel(y2)), y2, 'r-', 'LineWidth', lw2);
fill3(xgrid, zeros(1,numel(y2)), y2, 'k', 'EdgeColor', 'None', 'FaceColor', lgray);
axis([0,1,0,1,0,max(y3(:))]);
zlabel('Model evidence');
set(gca, 'Xgrid', 'Off', 'Ygrid', 'Off');

%% Customize the plots
%  ===================

for i = 1:9
    subplot(9,1,i);
    set(gca, 'XTick', [], 'YTick', [], 'ZTick', [], ...
        'TickDir', 'Out', 'TickLength', ones(1,2)./30, ...
        'LineWidth', lw2, 'FontSize', fontsize, 'FontWeight', 'Bold', ...
        'Box', 'Off', 'Layer', 'Bottom');
    set(get(gca, 'XLabel'), 'String', '');
    set(get(gca, 'YLabel'), 'String', '');
    set(get(gca, 'ZLabel'), 'String', '');
    xlim([0,1]);
end