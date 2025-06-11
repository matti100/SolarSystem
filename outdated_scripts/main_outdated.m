clear
clc
close all

%% CUSTOMIZE

% Flags
unit = 1;           % 1: day        Km/day
                    % 2: seconds    Km/s
         
frame = 1;          % 1: Solar system baricenter
                    % 0: Sun

plot_flag = 1;      % 1: plot
                    % 0: no plot

animation_flag = 0; % 1: animation
                    % 0: no animation

onlyEarth_flag = 0; % 1: plot only earth
                    % 0: plot all bodies

onlyIC_flag = 1;    % 1: extract only IC
                    % 0: extract all the ephemeris

sphere_flag = 0;    % 1: plot sphere
                    % 0: no sphere


% Time Interval for simulation
startTime = '2026-May-5';
stopTime = '2030-May-19';

startDateObj = datetime(startTime, 'InputFormat', 'yyyy-MMMM-d');
endDateObj = datetime(stopTime, 'InputFormat', 'yyyy-MMMM-d');

% Integration step
dt = 1; % [day]

% Plots
linewidth = 1.2;
markersize = 10;

%% PARAMETERS

% Inertial data
m_sun = 1.99 * 1e30;

m_mercury = 3.3 * 1e23;
m_venus = 4.87 * 1e24;
m_earth = 5.98 * 1e24;
m_mars = 6.42 * 1e23;
m_jupiter = 1.9 * 1e27;
m_saturn = 5.68 * 1e26;
m_uranus = 8.68 * 1e25;
m_neptune = 1.02 * 1e26;

M = [m_sun;
     m_mercury;
     m_venus;
     m_earth;
     m_mars;
     m_jupiter;
     m_saturn;
     m_uranus;
     m_neptune];

bodyNames = {'Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter',...
            'Saturn', 'Uranus', 'Neptune'};

bodyColors = {'#EDB120', '#D95319', '#77AC30', '#0072BD', '#A2142F',...
    '#7E2F8E', 'y', 'c', 'b'};

bodyRadius = [696000, 2440, 6052, 6371, 3390, 69911, 58232, 25362, 24622]; % [km]
scaleFactor = [20, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000];
caleFactor = ones(9,1);
bodySize = bodyRadius .* scaleFactor;
bodySize = [4, 1, 1, 1, 1, 2, 2, 2, 2];
[bodyR, bodyPhi, bodyTheta] = sphere(20);  % represent the single body (planet, sun ecc) as a sphere

% Universal gravitational constant
G_standard = 6.67430e-11; % [ m^3 / kg / s^2 ]
G = G_standard / (1000)^3 * (24*60*60)^2; % [ km^3 / kg / day^2 ]

% Astronomic Unit
AU = 149597870707 * 1e-3; % [km]

%% INITIALIZE PROBLEM

% Tspan
tspan = endDateObj - startDateObj;
t_len = days(tspan)/dt;

% Initial Condition (extract ephemeris)
[Ephemeris, x0, x_end, Trajectory, Velocity] = extractEphemeris(startTime, stopTime,...
                                                unit, dt, frame, onlyIC_flag);

% Initialize state vector 
N = length(M);  % Number of bodies
n = N*6;        % Number of integration variables [ (x, y, z, vx, vy, vz) * N ]
x = zeros(n, t_len);

%% SOLVE DIFFERENTIAL SYSTEM

% Update initial condition
x(:, 1) = x0;

% RK4 derivation method
%
%       x_i+1 = x_i + dt*sum_i(b_i * k_i)
%       ki = F(t + c_i*dt, x + dt*sum_j(a_ij * k_j)
%
%       with:   F = F(t, x)
%               a_ij, b_i, c_i from Butcher Tableau
%               
%       Butcher Tableau for RK4
%       | c1: 0                                                     |
%       | c2: 0.5      a21: 0.5                                     |
%       | c3: 0.5      a31: 0       a32: 0.5                        |
%       | c4: 1        a41: 0       a42: 0      a43: 1              |
%       |              b1: 1/6      b2: 1/3     b3: 1/3     b4: 1/6 |
% 
%       so:      
%       k1 = f( t, x )
%       k2 = f( t + 0.5*dt, x + dt*(0.5*k1) )
%       k3 = f( t + 0.5*dt, x + dt*(0.5*k2) )
%       k4 = f( t + dt, x + dt*(1*k3) )
% 
%       x_i+1 = x_i + dt*[ (1/6 * k1) + (1/3 * k2) + (1/3 * k3) + (1/6 * k4) ] 
%
% Note: in this case, the system is time-invariant 
%           f = f(x), not f = f(t, x)

% Apply Runge-Kutta 4
for i = 1:t_len-1
    k1 = BuildFunction(x(:, i), M, G);
    k2 = BuildFunction(x(:, i) + dt*0.5*k1, M, G);
    k3 = BuildFunction(x(:, i) + dt*0.5*k2, M, G);
    k4 = BuildFunction(x(:, i) + dt*1*k3, M, G);

    x(:, i+1) = x(:, i) + dt*( (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4);
end

%% PLOT
if plot_flag

    if onlyEarth_flag

        figure('Name', 'Plots - Earth');

        if onlyIC_flag == 0
            ax_ephem = subplot(2,1,1);
            hold(ax_ephem, 'on');
            % Trajectory
            plot3(Trajectory.Earth(1, :), Trajectory.Earth(2, :), Trajectory.Earth(3, :),...
                'LineWidth', linewidth, 'Color', bodyColors{4});
            plot3(Trajectory.Earth(1, end), Trajectory.Earth(2, end), Trajectory.Earth(3, end),...
                'Marker', 'o', 'MarkerEdgeColor', bodyColors{4}, 'MarkerFaceColor', bodyColors{4});

            if sphere_flag
            % Body
            % compute sphere
            [X, Y, Z] = createBody(bodyR, bodyPhi, bodyTheta, bodySize(4), ...
                [Trajectory.Earth(1, end), Trajectory.Earth(2, end), Trajectory.Earth(3, end)]);
            % plot sphere
            surf(X, Y, Z, 'FaceColor', bodyColors{4}, 'EdgeColor', 'none');
            % add name text
            text(Trajectory.Earth(1, end), Trajectory.Earth(2, end), Trajectory.Earth(3, end) + bodySize(4)*1.2, ...
                bodyNames{4}, 'Color', bodyColors{4}, ...
                'FontSize', 10, 'HorizontalAlignment', 'center');
            end
            hold(ax_ephem, 'off');
            view(ax_ephem, 3);

            % Options
            xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
            title('Trajectory - ephem', 'color', 'w');
            grid on;
            leg = legend('Earth trajectory');
            set(gca, 'Color', 'k');            % black background
            set(gcf, 'Color', 'k');            % black background of the figure
            set(gca, 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w'); % white axes
            set(gca, 'GridColor', 'w');        % white grid
            set(leg, 'TextColor', 'w');        % white legend
            set(leg, 'EdgeColor', 'w');
            set(leg, 'Color', 'k');
            if sphere_flag
            axis equal
            axis tight
            camlight
            lighting gouraud
            end
        end

        if onlyIC_flag == 0
            ax_sol = subplot(2,1,2);
        else
            ax_sol = axes('Position', [0.05 0.1 0.9 0.8]); % [left bottom width height]
        end
        idx = 6*(4-1)+1; % Earth index
        hold(ax_sol, 'on');
        % Trajectory
        plot3(x(idx, :), x(idx+1, :), x(idx+2, :), 'LineWidth', linewidth, ...
            'Color', bodyColors{4});
        plot3(x(idx, end), x(idx+1, end), x(idx+2, end), ...
            'Marker', 'o', 'MarkerEdgeColor', bodyColors{4}, 'MarkerFaceColor', bodyColors{4});
        
        if sphere_flag
        % Body
        % compute sphere
        [X, Y, Z] = createBody(bodyR, bodyPhi, bodyTheta, bodySize(4), ...
            [x(idx, end), x(idx+1, end), x(idx+2, end)]);
        % plot sphere
        surf(X, Y, Z, 'FaceColor', bodyColors{4}, 'EdgeColor', 'none');
        % add name text
        text(x(idx,end), x(idx+1,end), x(idx+2,end) + bodySize(4)*1.2, ...
            bodyNames{4}, 'Color', bodyColors{4}, ...
            'FontSize', 10, 'HorizontalAlignment', 'center');
        end
        hold(ax_sol, 'off');
        view(ax_sol, 3);

        % Options
        xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
        title('Trajectory - solved', 'color', 'w');
        grid on;
        leg = legend('Earth trajectory');
        set(gca, 'Color', 'k');            % black background
        set(gcf, 'Color', 'k');            % black background of the figure
        set(gca, 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w'); % white axes
        set(gca, 'GridColor', 'w');        % white grid
        set(leg, 'TextColor', 'w');        % white legend
        set(leg, 'EdgeColor', 'w');
        set(leg, 'Color', 'k');
        if sphere_flag
        axis equal
        axis tight
        camlight
        lighting gouraud
        end

    else

        figure('Name', 'Plots - all bodies');
        % axes('Position', [0.05 0.1 0.9 0.8]); % [left bottom width height]
        Traj = [];
        hold on;
        view(3);
        for i = 1:N
            j = 6*(i-1) + 1;
            % Trajectory
            traj = plot3(x(j, :), x(j+1, :), x(j+2, :), 'LineWidth', linewidth, ...
                'Color', bodyColors{i});
            plot3(x(j, end), x(j+1, end), x(j+2, end), ...
                'Marker', 'o', 'MarkerSize', bodySize(i), ...
                'MarkerEdgeColor', bodyColors{i}, 'MarkerFaceColor', bodyColors{i});
            Traj = [Traj; traj];
            
            if sphere_flag
            % Body
            [X, Y, Z] = createBody(bodyR, bodyPhi, bodyTheta, bodySize(i), ...
                [x(j, end), x(j+1, end), x(j+2, end)]);
            surf(X, Y, Z, 'FaceColor', bodyColors{i}, 'EdgeColor', 'none');
            text(x(j,end), x(j+1,end), x(j+2,end) + bodySize(i)*1.2, ...
                bodyNames{i}, 'Color', bodyColors{i}, ...
                'FontSize', 10, 'HorizontalAlignment', 'center');
            end
        end
        hold off;

        % Options
        title('Trajectory - all bodies');
        xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
        grid on;
        leg = legend(Traj, bodyNames);
        set(gca, 'Color', 'k');            % black background
        set(gcf, 'Color', 'k');            % black background of the figure
        set(gca, 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w'); % white axes
        set(gca, 'GridColor', 'w');        % white grid
        set(leg, 'TextColor', 'w');        % white legend
        set(leg, 'EdgeColor', 'w');
        set(leg, 'Color', 'k');
        if sphere_flag
        axis equal
        axis tight
        camlight
        lighting gouraud
        end
    end

end

%% ANIMATION
if animation_flag

    if onlyEarth_flag 

        idx = 6*(4-1)+1; % Earth index

        figure('Name', 'Animation - Earth');
        set(gcf, 'Color', 'k'); % black figure background
        
        if onlyIC_flag == 0
            % Ephemeris objects initializations
            ax_ephem = subplot(2,1,1);
            ax_sol = subplot(2,1,2);

            hold(ax_ephem, 'on');
            % Pre-allocate trajectory graphic object
            h_plot_ephem = plot3(ax_ephem, NaN, NaN, NaN, 'LineWidth', linewidth, 'Color', bodyColors{4});
            if sphere_flag
            % Pre-allocate Earth surface graphic object
            h_surf_ephem = surf(ax_ephem, NaN(2), NaN(2), NaN(2), 'FaceColor', bodyColors{4}, 'EdgeColor', 'none');
            % Pre-allocate text graphic object
            h_text_ephem = text(ax_ephem, NaN, NaN, NaN + bodySize(4)*1.2, ...
                            bodyNames{4}, 'Color', bodyColors{4}, ...
                            'FontSize', 10, 'HorizontalAlignment', 'center');
            end
            hold(ax_ephem, 'off');
            % Options - ephemeris
            xlabel(ax_ephem, 'X [km]'); ylabel(ax_ephem, 'Y [km]'); zlabel(ax_ephem, 'Z [km]');
            title(ax_ephem, 'Trajectory - ephem', 'color', 'w');
            stylePlot(ax_ephem, onlyEarth_flag, sphere_flag);
            % xmin = min(Trajectory.Earth(1, 1));
            % xmax = max(Trajectory.Earth(1, 1));
            % ymin = min(Trajectory.Earth(2, 1));
            % ymax = max(Trajectory.Earth(2, 1));
            % zmin = min(Trajectory.Earth(3, 1));
            % zmax = max(Trajectory.Earth(3, 1));
            % xlim(ax_ephem, [xmin xmax]);
            % ylim(ax_ephem, [ymin ymax]);
            % zlim(ax_ephem, [zmin zmax]);

            % Initialize trajectory vector;
            rx = [];
            ry = [];
            rz = [];

            view(ax_ephem, 3);
        else
            ax_sol = gca;
        end
        view(ax_sol, 3);

        % Solved solution objects initialization
        hold(ax_sol, 'on');
        % Pre-allocate trajectory graphic object
        h_plot_sol = plot3(ax_sol, NaN, NaN, NaN, 'LineWidth', linewidth, 'Color', bodyColors{4});
        if sphere_flag
        % Pre-allocate Earth surface graphic object
        h_surf_sol = surf(ax_sol, NaN(2), NaN(2), NaN(2), 'FaceColor', bodyColors{4}, 'EdgeColor', 'none');
        % Pre-allocate text graphic object
        h_text_sol = text(ax_sol, NaN, NaN, NaN + bodySize(4)*1.2, ...
                        bodyNames{4}, 'Color', bodyColors{4}, ...
                        'FontSize', 10, 'HorizontalAlignment', 'center');
        end
        hold(ax_sol, 'off');
        % Options - solved
        xlabel(ax_sol, 'X [km]'); ylabel(ax_sol, 'Y [km]'); zlabel(ax_sol, 'Z [km]');
        title(ax_sol, 'Trajectory - solved', 'color', 'w');
        stylePlot(ax_sol, onlyEarth_flag, sphere_flag);
        % xmin = min(x(idx, :));
        % xmax = max(x(idx, :));
        % ymin = min(x(idx +1, :));
        % ymax = max(x(idx +1, :));
        % zmin = min(x(idx +2, :));
        % zmax = max(x(idx +2, :));
        % xlim(ax_sol, [xmin xmax]);
        % ylim(ax_sol, [ymin ymax]);
        % zlim(ax_sol, [zmin zmax]);

        % Initialize trajectory vector
        xx = [];
        xy = [];
        xz = [];

        % Animate
        for i = 1:t_len
            
            % Animate Ephemeris
            if onlyIC_flag == 0
                
                % Update trajectory vector
                rx = [rx; Trajectory.Earth(1, i)];
                ry = [ry; Trajectory.Earth(2, i)];
                rz = [rz; Trajectory.Earth(3, i)];

                % Update body vector
                if sphere_flag
                [X, Y, Z] = createBody(bodyR, bodyPhi, bodyTheta, bodySize(4), ...
                    [rx(end), ry(end), rz(end)]);
                end

                % Update plot
                set(h_plot_ephem, 'XData', rx, 'YData', ry, 'ZData', rz);
                if sphere_flag
                set(h_surf_ephem, 'XData', X, 'YData', Y, 'ZData', Z);
                set(h_text_ephem, 'Position', [rx(end), ry(end), rz(end) + bodySize(4)*1.2]);
                end
            end
           
            % Animate solved solution
            % Update trajectory vector
            xx = [xx; x(idx, i)];
            xy = [xy; x(idx+1, i)];
            xz = [xz; x(idx+2, i)];

            % Update body vector
            [X, Y, Z] = createBody(bodyR, bodyPhi, bodyTheta, bodySize(4), ...
                [xx(end), xy(end), xz(end)]);
            
            % Update plot
            set(h_plot_sol, 'XData', xx, 'YData', xy, 'ZData', xz);
            if sphere_flag
            set(h_surf_sol, 'XData', X, 'YData', Y, 'ZData', Z);
            set(h_text_sol, 'Position', [xx(end), xy(end), xz(end) + bodySize(4)*1.2]);
            end

            drawnow;

            % Pause for a better visualization
            pause(0.005);
        end

    else
        %% 
    end
end


%% Functions
function [X, Y, Z] = createBody(r, phi, theta, size, pos)

X = size * r + pos(1);
Y = size * phi + pos(2);
Z = size * theta + pos(3);

end

function stylePlot(ax, onlyEarth_flag, sphere_flag)

% Options for plot
grid on;
if onlyEarth_flag
    leg = legend(ax, 'Earth trajectory');
else 
    leg = legend(ax, bodyNames);
end
set(ax, 'Color', 'k');            % black background
set(ax, 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w'); % white axes
set(ax, 'GridColor', 'w');        % white grid
set(leg, 'TextColor', 'w');        % white legend
set(leg, 'EdgeColor', 'w');
set(leg, 'Color', 'k');

if sphere_flag
camlight(ax)
lighting(ax, 'gouraud');
axis(ax, 'equal');
axis(ax, 'tight');
end

end