% This script simulates choices of an observer can take in the the
% Monty-Hall problem. There are two possible policies to adopt once the TV
% presenter opened a door behind which there is a goat: one can either (i)
% stick to his/her first choice  or (ii) change to (one of) the remaining
% door(s). This simulation shows (as theory) that the second policy is
% actually better. This is true no matter of the total number of doors in
% the game.
%
% Maxime Maheu, 11/2016

%% Initialization
%  ==============

% Clear the place
clear; close('all'); clc;

% Set some options
iterplot = true;
dispoutp = true;

% Define the number of problems to simulate
Nsimu = 100;

% Define the number of doors (usually it is set to 3)
Ndoors = 3;
doors = 1:Ndoors;

% Initialize random generators
rng('shuffle');

% Prepare outputs
DumbWins  = NaN(1,Nsimu);
SmartWins = NaN(1,Nsimu);

% Prepare the figure
F1 = figure;
x = 0:0.001:1;
set(gca, 'XTick', 0:0.1:1, 'LineWidth', 1, 'FontSize', 18, 'Box', 'Off', 'Layer', 'Bottom');
xlabel('Probability of winning the car');
ylabel('Probability mass');

% Define possible outcomes
Outcomes = {'Fail', 'Win'};

%% Run the loop
%  ============

% For each simulation
for i = 1:Nsimu
    
    %% Set-up the prices behind the doors
    %  ==================================
    
    % Hide teh car behind one of the door
    CarDoor = randi(Ndoors);
    
    % Hide goats behind the remaining doors
    GoatDoors = setdiff(doors, CarDoor);
    
    %% First round: pick a door
    %  ========================
    
    % Make both observer choose the same door
    FirstChoosenDoor = randi(Ndoors);
    
    % Make the TV presenter opens one of the door that is hidding a goat
    RevealedDoor = GoatDoors(randi(Ndoors-1));
    
    %% Second round: change your choice or not
    %  =======================================
    
    % The first observer sticks to his/her first choice
    SecondChoosenDoor_Dumb = FirstChoosenDoor;
    
    % The second observer changes his/her choice to one of the remaining
    % door (i.e. one that is not open or previously choosen)
    RemainingDoors = setdiff(doors, [FirstChoosenDoor, RevealedDoor]);
    SecondChoosenDoor_Smart = RemainingDoors(randi(numel(RemainingDoors)));
    
    %% Reveal the door that was hiding the car and see whether observers have won or not
    %  =================================================================================
    
    % Did the dumb player chose the door hiding the car?
    DumbWins(i)  = SecondChoosenDoor_Dumb  == CarDoor;
    
    % Did the smart player chose the door hiding the car?
    SmartWins(i) = SecondChoosenDoor_Smart == CarDoor;
    
    %% Display the simulation in a figure
    %  ==================================
    
    if iterplot || i == Nsimu
        
        % Derive beliefs in the probability of winning the car
        D = betapdf(x, nansum(DumbWins  == 1) + 1, nansum(DumbWins  == 0) + 1);
        S = betapdf(x, nansum(SmartWins == 1) + 1, nansum(SmartWins == 0) + 1);
        
        % Normalize it
        D = D ./ nansum(D);
        S = S ./ nansum(S);
        
        % Display the results
        F1; cla; hold('on');
        fill([0,x,1], [0,D,0], 'b-', 'LineWidth', 2, 'FaceAlpha', 0.5);
        fill([0,x,1], [0,S,0], 'r-', 'LineWidth', 2, 'FaceAlpha', 0.5);
        plot(1/Ndoors, 0, 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 15, 'LineWidth', 2);
        plot(1/(Ndoors-1), 0, 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 15, 'LineWidth', 2);
        legend({'Dumb player: always stick to its first choice', ...
                'Smart player: always shift to (one of) the remaining door(s)'}, ...
                'Location', 'SouthOutside');
        xlim([0,1]);
        drawnow;
    end
    
    %% Display the output in the command window
    %  ========================================
    
    if dispoutp
        fprintf('- Simulation %3i: ', i);
        fprintf('dumb player %ss (%g%%), ', lower(Outcomes{DumbWins(i)+1}), nanmean(DumbWins  == 1)*100);
        fprintf('smart player %ss (%g%%).\n', lower(Outcomes{SmartWins(i)+1}), nanmean(SmartWins  == 1)*100);
    end
end