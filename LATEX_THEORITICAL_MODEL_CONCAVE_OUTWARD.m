% V SHAPED DOUBLE PEELING ON CONCAVE SUBSTRATE

clear; clc; close all;

% Given constants
R = 5;                    %  
mu = 0.1;                 % shear modulus in MPa
theta = pi/3;             
t = 0.1;                  % thicknes is 0.1 mm  

% Define actual G (shear modulus of adhesive or layer)
G_values_ini = [0, 0.0005, 0.005, 0.025, 0.05];  % in MPa or consistent units

% Compute corresponding dimensionless Gbar
Gbar_values = 2 * G_values_ini / (mu * t);

% Define Gbar values
Gbar_values = [0, 0.1, 1, 5, 10];    % defining various Gbar values we want to analyse 

% Range of ubar values
ubar_values = linspace(0, 5, 5000);

% Solver options
options = optimset('Display', 'off', 'TolFun', 1e-10, 'TolX', 1e-10);      % suppressing iterative output & setting up the function value tolerance & solution tolerance 

% Preallocate a matrix to store solutions
L1bar_solutions = zeros(length(Gbar_values), length(ubar_values));
convergence_flags = zeros(length(Gbar_values), length(ubar_values));

% Function F(Gbar, L1bar, ubar) with better numerical handling
F = @(Gbar, L1bar, ubar) ...                                               % anonymous function F with three inputs: (Gbar, L1bar, ubar)
    (Gbar + (2.*(1-ubar).*sin(L1bar)./(L1bar + 1e-12)) ...
    + (4*L1bar ./ sqrt(1e-12 + (-1 + (ubar + cos(L1bar))).^2 + (sin(L1bar)).^2)) ...
    + ((-2 + 2.*ubar - ubar.^2 + 2.*((1-ubar).*cos(L1bar)))./(L1bar.^2 + 1e-12)) ...
    - 3 ...
    - (2.*L1bar.^2.*(1-ubar).*sin(L1bar) ./ ((1e-12 + ((-1 + (ubar + cos(L1bar))).^2 + (sin(L1bar)).^2)).^(3/2))));

%% SOLVER SECTION

numG = numel(Gbar_values);
numU = numel(ubar_values);

L1bar_solutions = nan(numG,numU);
conv_flags = zeros(numG,numU);

initial_guesses = ...
[1e-12 1e-10 1e-8 1e-6 1e-4 1e-3 1e-2 0.05 0.1];

for ig = 1:numG

    Gbar = Gbar_values(ig);

    fprintf('\nSolving for Gbar = %.2f\n',Gbar);

    for iu = 1:numU

        u = ubar_values(iu);

        F_L1 = @(L1) F(Gbar,L1,u);

        roots_found = [];

        %% Continuation guess

        if iu == 1

            guess_list = initial_guesses;

        else

            prev = L1bar_solutions(ig,iu-1);

            if ~isnan(prev)

                guess_list = [ ...
                    prev ...
                    0.5*prev ...
                    2*prev ...
                    initial_guesses ];

            else

                guess_list = initial_guesses;

            end

        end

        %% Search using fsolve

        for guess = guess_list

            try

                [L1sol,fval,exitflag] = ...
                    fsolve(F_L1,guess,options);

                if exitflag > 0 && ...
                   abs(fval) < 1e-8 && ...
                   isreal(L1sol) && ...
                   L1sol > 0 && ...
                   L1sol < 10

                    roots_found(end+1) = L1sol;

                end

            catch
            end

        end

        %% Keep smallest positive root

        if ~isempty(roots_found)

            roots_found = unique(round(roots_found,10));

            L1bar_solutions(ig,iu) = min(roots_found);

            conv_flags(ig,iu) = 1;

        else

            %% Sign-change scan fallback

            L1scan = logspace(-12,1,10000);

            Fscan = arrayfun(F_L1,L1scan);

            idx = find(diff(sign(Fscan)));

            if ~isempty(idx)

                guess2 = L1scan(idx(1));

                try

                    [L1sol,fval,exitflag] = ...
                        fsolve(F_L1,guess2,options);

                    if exitflag > 0 && ...
                       abs(fval) < 1e-8 && ...
                       L1sol > 0

                        L1bar_solutions(ig,iu) = L1sol;

                        conv_flags(ig,iu) = 2;

                    end

                catch
                end

            end

        end

    end

end

%% Plot L1bar vs ubar

figure; hold on; grid on; box on;
cols = lines(numG);

L1bar_solutions(:,1) = 0;
conv_flags(:,1) = 1;

for ig = 1:numG
    ok = conv_flags(ig,:)>0;
    
    % Lplot = [0, L1bar_solutions(ig,ok)];
    % uplot = [0, ubar_values(ok)];

    plot(ubar_values(ok), L1bar_solutions(ig,ok), 'Color',cols(ig,:), 'LineWidth',4, ...
        'DisplayName',sprintf('\\bar{G}=%.2f',Gbar_values(ig)));
end
xlabel('$\bar{u}$','Interpreter','latex');
ylabel('$\bar{L}_1$','Interpreter','latex');
% title('$\bar{L}_1$ vs $\bar{u}$','Interpreter','latex');
legend('Location','best','Interpreter','latex');

set(gca,'FontSize',35)
xlabel('Dimensionless vertical displacement, $\bar{u}$','Interpreter','latex','FontSize', 35);
ylabel('Dimensionless peeling length, $\bar{L_1}$','Interpreter','latex','FontSize', 35);
% title('$\bar{L}_1$ vs $\bar{u}$','Interpreter','latex','FontSize', 30);
legend('Location','best','Interpreter','latex');
lgd = legend({'$\bar{G}=0.0$','$\bar{G}=0.1$','$\bar{G}=1.0$','$\bar{G}=5.0$','$\bar{G}=10.0$'}, ...
             'Interpreter','latex');
legend boxoff
legend('FontSize',35)
ylim([0 1.6]);
xlim([0 5]);
ax = gca;
ax.GridLineStyle = '--';
ax.GridColor = [0 0 0];

hold on

%% Compute ALPHA
alpha_solutions = nan(size(L1bar_solutions));
for ig = 1:numG
    for iu = 1:numU
        L1 = L1bar_solutions(ig, iu);
        if ~isnan(L1)
            u = ubar_values(iu);
            alpha_solutions(ig, iu) = pi/2 + L1 - atan2(sin(L1), (u - 1 + cos(L1)));   % FORMULA FOR ALPHA
        end
    end
end

% Plot alpha vs ubar
figure; hold on; grid on; box on;
for ig = 1:numG
    ok = conv_flags(ig,:)>0;
    
    % Assign specific colors based on exact Gbar value
    switch Gbar_values(ig)
        case 0
            line_color = [0.7 0.7 0.7];  % grey
        case 0.1
            line_color = [1 0 0];        % red
        case 1
            line_color = [0 0 1];        % blue
        case 5
            line_color = [0 1 0];        % green
        case 10
            line_color = [0.5 0 0.5];    % purple
        otherwise
            line_color = [0 0 0];        % black for others
    end
    
    plot(ubar_values(ok), alpha_solutions(ig,ok), 'Color', line_color, 'LineWidth', 4, ...
        'DisplayName', sprintf('\\bar{G}=%.1f', Gbar_values(ig)));
end
% xlabel('$\bar{u}$','Interpreter','latex');             % 'Dimensionless vertical displacement, $\bar{u}$', 'Interpreter', 'latex'
% ylabel('$\alpha$','Interpreter','latex');
% title('$\alpha$ vs $\bar{u}$','Interpreter','latex');
% legend('Location','best','Interpreter','latex');
% ylim([0 3.5]);

set(gca,'FontSize',35)
xlabel('Dimensionless vertical displacement, $\bar{u}$','Interpreter','latex','FontSize', 35);
ylabel('Peeling angle, $\alpha$','Interpreter','latex','FontSize', 35);
% title('$\bar{L}_1$ vs $\bar{u}$','Interpreter','latex','FontSize', 30);
legend('Location','best','Interpreter','latex');
lgd = legend({'$\bar{G}=0.0$','$\bar{G}=0.1$','$\bar{G}=1.0$','$\bar{G}=5.0$','$\bar{G}=10.0$'}, ...
             'Interpreter','latex');
legend boxoff
legend('FontSize',35)
ylim([0 3]);
xlim([0 5]);
ax = gca;
ax.GridLineStyle = '--';
ax.GridColor = [0 0 0];


% % marking red coloured reference lines
% yline(3, '--r', 'L_1 = 3', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
% xline(0.2525, '--r', 'L_1 = 3', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');




%% Compute BETA

beta_solutions = nan(size(L1bar_solutions));
for ig = 1:numG
    for iu = 1:numU
        L1 = L1bar_solutions(ig, iu);
        if ~isnan(L1)
            u = ubar_values(iu);
            beta_solutions(ig, iu) = atan2(sin(L1), (u - 1 + cos(L1)));    % FORMULA FOR BETA 
        end
    end
end

% Plot BETA vs ubar
figure; hold on; grid on; box on;
for ig = 1:numG
    ok = conv_flags(ig,:)>0;
    
    % Assign specific colors based on exact Gbar value
    switch Gbar_values(ig)
        case 0
            line_color = [0.7 0.7 0.7];  % grey
        case 0.1
            line_color = [1 0 0];        % red
        case 1
            line_color = [0 0 1];        % blue
        case 5
            line_color = [0 1 0];        % green
        case 10
            line_color = [0.5 0 0.5];    % purple
        otherwise
            line_color = [0 0 0];        % black for others
    end
    
    plot(ubar_values(ok), beta_solutions(ig,ok), 'Color', line_color, 'LineWidth', 4, ...
        'DisplayName', sprintf('\\bar{G}=%.1f', Gbar_values(ig)));
end

% xlabel('$\bar{u}$','Interpreter','latex');
% ylabel('$\beta (rad)$','Interpreter','latex');
% title('$\beta (rad)$ vs $\bar{u}$','Interpreter','latex');
% legend('Location','best','Interpreter','latex');
% ylim([0 3.5]);

set(gca,'FontSize',35)
xlabel('Dimensionless vertical displacement, $\bar{u}$','Interpreter','latex','FontSize', 35);
ylabel({'Angle between vertical peeling force'; ...
       'and peeling arm, $\beta$'}, ...
      'Interpreter','latex', ...
      'FontSize',35);
% title('$\bar{L}_1$ vs $\bar{u}$','Interpreter','latex','FontSize', 30);
legend('Location','best','Interpreter','latex');
lgd = legend({'$\bar{G}=0.0$','$\bar{G}=0.1$','$\bar{G}=1.0$','$\bar{G}=5.0$','$\bar{G}=10.0$'}, ...
             'Interpreter','latex');
legend boxoff
legend('FontSize',35)
ylim([0 1.6]);
xlim([0 5]);
ax = gca;
ax.GridLineStyle = '--';
ax.GridColor = [0 0 0];

hold on

% yline(3, '--r', 'L_1 = 3', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
% xline(0.2525, '--r', 'L_1 = 3', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
% yline(0.809079, '--r', 'L_1 = 3', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
% 


%% Compute STRETCH - CORRECTED VERSION

for ig = 1:numG
    fprintf('Gbar = %.1f, min(L1) = %.3e\n', ...
        Gbar_values(ig), ...
        min(L1bar_solutions(ig,:),[],'omitnan'));
end


stretch_solutions = nan(size(L1bar_solutions));

for ig = 1:numG
    for iu = 1:numU

        L1 = L1bar_solutions(ig, iu);

        if ~isnan(L1)

            u = ubar_values(iu);

            if abs(L1) < 1e-14
                stretch_solutions(ig, iu) = 1;
            else
                stretch_solutions(ig, iu) = ...
                    sqrt((u - 1 + cos(L1))^2 + sin(L1)^2)/L1;
            end

        end 

    end

end

% Plot STRETCH vs ubar
figure; hold on; grid on; box on;
for ig = 1:numG
    ok = conv_flags(ig,:)>0;
    
    % Assign specific colors based on exact Gbar value
    switch Gbar_values(ig)
        case 0
            line_color = 'k';  % grey
        case 0.1
            line_color = [0.7 0.7 0.7];        % red
        case 1
            line_color = [1 0 0];        % blue
        case 5
            line_color = [0 0 1];        % green
        case 10
            line_color = [0 1 0];        % purple
        otherwise
            line_color = 'k';        % black for others
    end
    
    plot(ubar_values(ok), stretch_solutions(ig,ok), 'Color', line_color, 'LineWidth', 4, ...
        'DisplayName', sprintf('\\bar{G}=%.1f', Gbar_values(ig)));
end

% xlabel('$\bar{u}$','Interpreter','latex');
% ylabel('$\lambda$','Interpreter','latex');
% title('$\lambda$ vs $\bar{u}$','Interpreter','latex');
% legend('Location','best','Interpreter','latex');
% ylim([0 3.5]);

set(gca,'FontSize',35)
xlabel('Dimensionless vertical displacement, $\bar{u}$','Interpreter','latex','FontSize', 35);
ylabel('Stretch of peeling arm, $\lambda$','Interpreter','latex','FontSize', 35);
% title('$\bar{L}_1$ vs $\bar{u}$','Interpreter','latex','FontSize', 30);
legend('Location','best','Interpreter','latex');
lgd = legend({'$\bar{G}=0.0$','$\bar{G}=0.1$','$\bar{G}=1.0$','$\bar{G}=5.0$','$\bar{G}=10.0$'}, ...
             'Interpreter','latex');
legend boxoff
legend('FontSize',35)
ylim([1 3.5]);
xlim([0 5]);
ax = gca;
ax.GridLineStyle = '--';
ax.GridColor = [0 0 0];


hold on


%% Compute TENSION
Tension_solutions = nan(size(L1bar_solutions));
for ig = 1:numG
    for iu = 1:numU
        L1 = L1bar_solutions(ig, iu);
        if ~isnan(L1)
            u = ubar_values(iu);
            % CORRECTED TENSION FORMULA (dimensionless)
            Tension_solutions(ig, iu) = (stretch_solutions(ig, iu) - (1 ./ (stretch_solutions(ig, iu).^2)));  % Note : By multiplying R to sin(L1) the plot is shifted to topside.
        end 
    end
end

% Plot TENSION vs ubar
figure; hold on; grid on; box on;
for ig = 1:numG
    ok = conv_flags(ig,:)>0;
    
    % Assign specific colors based on exact Gbar value
    switch Gbar_values(ig)
        case 0
            line_color = 'k';  % grey
        case 0.1
            line_color = [0.7 0.7 0.7];        % red
        case 1
            line_color = [1 0 0];        % blue
        case 5
            line_color = [0 0 1];        % green
        case 10
            line_color = [0 1 0];        % purple
        otherwise
            line_color = 'k';        % black for others
    end
    
    plot(ubar_values(ok), Tension_solutions(ig,ok), 'Color', line_color, 'LineWidth', 4, ...
        'DisplayName', sprintf('\\bar{G}=%.1f', Gbar_values(ig)));
end

% xlabel('$\bar{u}$','Interpreter','latex');
% ylabel('$\bar{T}$','Interpreter','latex');
% title('$\bar{T}$ vs $\bar{u}$','Interpreter','latex');
% legend('Location','best','Interpreter','latex');
% ylim([0 3.5]);

% for Latex visuals
set(gca,'FontSize',35)
xlabel('Dimensionless vertical displacement, $\bar{u}$','Interpreter','latex','FontSize', 35);
ylabel('Dimensionless vertical peeling force, $\bar{T}$','Interpreter','latex','FontSize', 35);
% title('$\bar{L}_1$ vs $\bar{u}$','Interpreter','latex','FontSize', 30);
legend('Location','best','Interpreter','latex');
lgd = legend({'$\bar{G}=0.0$','$\bar{G}=0.1$','$\bar{G}=1.0$','$\bar{G}=5.0$','$\bar{G}=10.0$'}, ...
             'Interpreter','latex');
legend boxoff
legend('FontSize',35)
ylim([0 3.5]);
xlim([0 5]);
ax = gca;
ax.GridLineStyle = '--';
ax.GridColor = [0 0 0];
hold on

%% Compute FORCE - (DIMENSIONLESS) 
Force_solutions = nan(size(L1bar_solutions));
for ig = 1:numG
    for iu = 1:numU
        L1 = L1bar_solutions(ig, iu);
        if ~isnan(L1)
            u = ubar_values(iu);
            % CORRECTED FORCE FORMULA (dimensionless)
            Force_solutions(ig, iu) = cos(beta_solutions(ig, iu)).* Tension_solutions(ig, iu);   % Note : By multiplying R to sin(L1) the plot is shifted to topside.
        end 
    end
end

% Plot FORCE vs ubar
figure; hold on; box on; grid on;
for ig = 1:numG
    ok = conv_flags(ig,:)>0;
    
    % Assign specific colors based on exact Gbar value
    switch Gbar_values(ig)
        case 0
            line_color = 'k';  % grey
        case 0.1
            line_color = [0.7 0.7 0.7];        % red
        case 1
            line_color = [1 0 0];        % blue
        case 5
            line_color = [0 0 1];        % green
        case 10
            line_color = [0 1 0];        % purple
        otherwise
            line_color = 'k';        % black for others
    end
    
    plot(ubar_values(ok), Force_solutions(ig,ok), 'Color', line_color, 'LineWidth', 4, ...
        'DisplayName', sprintf('\\bar{G}=%.1f', Gbar_values(ig)));
end

% for Latex visuals 
set(gca,'FontSize',35)
xlabel('Dimensionless vertical displacement, $\bar{u}$','Interpreter','latex','FontSize', 35);
ylabel('Dimensionless vertical peeling force, $\bar{F}$','Interpreter','latex','FontSize', 35);
% title('$\bar{L}_1$ vs $\bar{u}$','Interpreter','latex','FontSize', 30);
legend('Location','best','Interpreter','latex');
lgd = legend({'$\bar{G}=0.0$','$\bar{G}=0.1$','$\bar{G}=1.0$','$\bar{G}=5.0$','$\bar{G}=10.0$'}, ...
             'Interpreter','latex');
legend boxoff
legend('FontSize',35)
ylim([0 3.5]);
xlim([0 5]);
ax = gca;
ax.GridLineStyle = '--';
ax.GridColor = [0 0 0];

hold on


%==================================================================================================================
% END OF CONCAVE FORMULATION LAST UPDATED 6 JULY 2026
%==================================================================================================================