% Toy example about how to compare general linear models using variational
% Bayesian inference as implemented in the VBA toolbox. In particular we
% show (i) how this approach handles models with different number of
% parameters and (ii) how random-effect is taken into account in group
% studies.
%
% Maxime Maheu, 06/2018

%% Simulate data from a group of subjects
%  ======================================

% For the sake of reproducbility, initialize the random generator
rng('default');

% Define dimensionality of the problem
nObs   = 200; % data dimension
nParam = [1, 3]; % parameters dimension
SubMod = [ones(1,17), ones(1,3)+1]; % repartition of models across subjects
nSub   = numel(SubMod); % number of subjects

% Define SNR options
nv = 1; % noise variance
sp = 1; % signal power

% Prepare the output variables
X    = NaN(nObs,max(nParam)+1,nSub);
Yhat = NaN(nObs,nSub,2);
Y    = NaN(nObs,nSub);

% For each subject
for iSub = 1:nSub
    
    % Create ubject-specific design matrix
    X(:,:,iSub) = [ones(nObs, 1), randn(nObs, max(nParam))];
    
    % Define subject-specific beta parameters
    b = sqrt(sp) * rand(max(nParam)+1, 1);
    
    % Simulate observations based on the design matrix and the (hidden)
    % beta parameters
    for iMod = 1:2
        Yhat(:,iSub,iMod) = X(:,1:nParam(iMod)+1,iSub) * b(1:nParam(iMod)+1);
    end

    % Generate some noise
    noise = sqrt(nv) * randn(nObs,1);
    
    % Generate simulated observations based on all or serveral predictors
    % depending on the subject (i.e. depending on the chosen generative
    % model for that particular subject)
    Y(:,iSub) = Yhat(:,iSub,SubMod(iSub)) + noise;
end

% Display simulated data
figure('Name', 'Simulated data', 'Units', 'Normalized', 'Position', [0.2 0.2 0.3 0.6]);
for iSub = 1:nSub
    
    % Data versus preditions from model 1
    subplot(4,5,iSub);
    yyaxis('left');
    plot(Y(:,iSub), Yhat(:,iSub,1), '.'); lsline;
    ylim([min(Yhat(:)), max(Yhat(:))]);
    
    % Data versus preditions from model 2
    yyaxis('right');
    plot(Y(:,iSub), Yhat(:,iSub,2), '.'); lsline;
    ylim([min(Yhat(:)), max(Yhat(:))]);
    xlim([min(Y(:)), max(Y(:))]);
end

%% Fit two simple models to the data of each subject
%  =================================================

% Whether to infer the priors from the data:
%   - >> empbayes = true; % estimate priors using empirical Bayes
%   - >> empbayes = false; % rely on user-defined priors
%   - >> empbayes = []; % rely on the frequentist limit 
empbayes = [];

% Rely on the frequentist limit to estimate the log-evidence of a GLM
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if isempty(empbayes)
    
% Prepare the output variable
L = NaN(2,nSub);

% For each subject
for iSub = 1:nSub
    
    % For each model
    for iMod = 1:2
        
        % Get the number of parameters of that model
        np = nParam(iMod) + 1;
        
        % Get (simulated) observations from that subject
        y = Y(:,iSub);
        
        % Get (simulated) subject's design matrix
        x = X(:,1:nParam(iMod)+1,iSub);
        
        % Use a frequentist limit to estimate the log-evidence
        L(iMod,iSub) = lev_GLM(y, x);
    end
end

% Use manually defined priors
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif ~empbayes

% Prepare output variables
post = cell(2,nSub); % posterior over parameters
out  = cell(2,nSub); % quality of fit

% For each subject
for iSub = 1:nSub
    
    % Observations to be explained
    y = Y(:,iSub);
    
    % For each model
    for iMod = 1:2
        np = nParam(iMod)+1;
        
        % Define predictions 
        options                 = [];
        options.inG.X           = X(:,1:np,iSub);
        
        % Whether to display 
        options.DisplayWin      = 0;
        options.verbose         = 1;
        
        % Define prior beliefs
        options.priors.muPhi    = ones(np,1)./2;
        options.priors.SigmaPhi = eye(np);
        
        % Define the number of evolution/observation parameters to estimate
        dim                     = [];
        dim.n_phi               = np;
        dim.n_theta             = 0;
        dim.n                   = 0;
        
        % Run the fitting routine
        [post{iMod,iSub}, out{iMod,iSub}] = ...
            VBA_NLStateSpaceModel(y, [], [], @g_GLM, dim, options);
    end
end

% Use free energy as a 
L = cellfun(@(x) x.F, out);

% Estimate priors using empirical Bayes
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif empbayes

% Data to explain
y = mat2cell(Y, nObs, ones(1, nSub));

% Prepare output variables
p_sub = cell(1,2); p_gp = cell(1,2); % posterior over parameters
o_sub = cell(1,2); o_gp = cell(1,2); % quality of fit

% For each model
for iMod = 1:2
    np = nParam(iMod)+1;
        
    % Specify the dimensions of the problem
    dim         = []; % empty variable
    dim.n_phi   = np; % number of observation parameters
    dim.n       = 0;  % number of hidden states
    dim.n_theta = 0;  % number of evolution parameters
    dim.n_t     = 1;  % number of time samples
    
    % Define minimal options
    options             = [];
    options.DisplayWin  = 0;
    options.verbose     = 0;
    options             = repmat({options}, [nSub,1]);
    optiongp            = [];
    optiongp.DisplayWin = 0;
    optiongp.verbose    = 1;
    
    % Explaining variable: position of the observation post-change-point
    for iSub = 1:nSub, options{iSub}.inG.X = X(:,1:np,iSub); end
    
    % Run the fitting routine while infering priors from the data using
    % empirical Bayes
    [p_sub{iMod}, o_sub{iMod}, p_gp{iMod}, o_gp{iMod}] = ...
        VBA_MFX(y, [], [], @g_GLM, dim, options, [], optiongp);
end

% Use the free energy as an approximation of the log-evidence
out = [o_sub{:}]';
L = cellfun(@(x) x.F, out);

end

%% Perform a Bayesian Model Selection to identify the best model across subjects
%  =============================================================================

% Perform the BMS using random effects
options = [];
options.DisplayWin = 0;
[posterior, output] = VBA_groupBMC(L, options);

% Compute the protected exceedance probabilities
output.pep = (1-output.bor)*output.ep + output.bor/length(output.ep);

%% Compare RFX and FFX bayesian inference
%  ======================================

% Prepare a new window
figure('Name', 'BMS results', 'Units', 'Normalized', 'Position', [0.5 0.2 0.3 0.6]);

% Random effects
% ~~~~~~~~~~~~~~

% Subject-specific model attribution
subplot(3,2,1);
barh(1:nSub, posterior.r', 1, 'stacked');
axis([0, 1, 1/2, nSub+1/2]);
xlabel('p(M_i|y)');
ylabel('Subject #');
legend({'M_1', 'M_2'});
title({'RFX', 'Model attribution'});

% Dirichlet distribution
subplot(3,2,3);
pgrid = linspace(0, 1, 1001);
dist = betapdf(pgrid, posterior.a(1), posterior.a(2));
fill([pgrid(pgrid >= 1/2), fliplr(pgrid(pgrid >= 1/2))], ...
    [dist(pgrid >= 1/2), zeros(1, sum(pgrid >= 1/2))], 'k', ...
    'FaceColor', ones(1,3)./2); hold('on');
plot(pgrid, dist, 'k-');
plot(ones(1,2)/2, ylim, 'k--');
set(gca, 'Box', 'Off');
xlabel('r_1');
ylabel('p(r_1|y)');
title('Dirichlet density describing the probability of M_1');

% Model frequencies and (protected) exceedance probabilities
subplot(3,2,5);
m = mean(posterior.r, 2);
nv = std(posterior.r, [], 2) ./ sqrt(nSub);
yyaxis('left');
bar(1:2, m); hold('on');
plot(repmat(1:2, 2, 1), (m+nv.*[-1,1])', 'k-');
plot([0,3], ones(1,2)./2, '--');
set(gca, 'XTick', 1:2, 'XTickLabel', {'M_1', 'M_2'}, 'Box', 'Off');
axis([0,3,0,1]);
ylabel('Model frequencies');
yyaxis('right'); lgd = NaN(1,2);
plot([0,3], repmat(0.95,1,2), '--'); hold('on');
lgd(1) = plot((1:2)-1/5, output.ep,  'o', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
lgd(2) = plot((1:2)+1/5, output.pep, '^', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
axis([0,3,0,1]);
legend(lgd, {'EP', 'PEP'});
ylabel('(protected) exceedance probabilities');
title('Key quantities');

% Fixed effects
% ~~~~~~~~~~~~~

% Individual bayes factors
subplot(3,2,2);
logBF = L(1,:)-L(2,:);
barh(1:nSub, logBF, 1, 'FaceColor', ones(1,3)./2);
set(gca, 'Box', 'Off');
axis([-max(abs(xlim)), max(abs(xlim)), 1/2, nSub+1/2]);
xlabel('log p(y|M_1) - log p(y|M_2)');
ylabel('Subject #');
title({'FFX', 'Individual (log) Bayes factor'});

% Group bayes factor
subplot(3,2,4);
m = mean(logBF);
nv = std(logBF) / sqrt(nSub);
bar(1, m, 'k'); hold('on');
plot(ones(1,2), m+nv*[-1,1], 'Color', ones(1,3)./2);
set(gca, 'XTick', [], 'Box', 'Off');
ylabel('Group-average log Bayes factor');
title('Global (log) Bayes factor');
