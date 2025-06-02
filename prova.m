clc; clear; close all;

% Costanti fisiche
G = 6.67430e-11; % Costante gravitazionale [m^3 kg^-1 s^-2]
M = 1.989e30;    % Massa del Sole [kg]

% Condizioni iniziali (posizione e velocità della Terra)
AU = 1.496e11;   % Unità astronomica [m] (distanza media Terra-Sole)
r0 = [AU; 0];    % Posizione iniziale [m]
v0 = [0; 29780]; % Velocità iniziale [m/s] (circa 29.78 km/s)

% Parametri di simulazione
dt = 60 * 60 * 24;  % Passo temporale [s] (1 giorno)
T = 365 * dt;       % Tempo totale di simulazione (1 anno)
N = floor(T / dt);  % Numero di passi

% Preallocazione array
r = zeros(2, N);
v = zeros(2, N);
t = zeros(1, N);

% Condizioni iniziali
r(:, 1) = r0;
v(:, 1) = v0;

% Funzione per accelerazione gravitazionale
acc = @(r) -G * M / norm(r)^3 * r;

% Simulazione con metodo di Runge-Kutta 4° ordine
for i = 1:N-1
    k1v = dt * acc(r(:, i));
    k1r = dt * v(:, i);

    k2v = dt * acc(r(:, i) + k1r / 2);
    k2r = dt * (v(:, i) + k1v / 2);

    k3v = dt * acc(r(:, i) + k2r / 2);
    k3r = dt * (v(:, i) + k2v / 2);

    k4v = dt * acc(r(:, i) + k3r);
    k4r = dt * (v(:, i) + k3v);

    v(:, i+1) = v(:, i) + (k1v + 2*k2v + 2*k3v + k4v) / 6;
    r(:, i+1) = r(:, i) + (k1r + 2*k2r + 2*k3r + k4r) / 6;
    
    t(i+1) = t(i) + dt;
end

% Grafico dell'orbita
figure;
plot(r(1, :), r(2, :), 'b', 'LineWidth', 1.5); hold on;
plot(0, 0, 'yo', 'MarkerSize', 12, 'MarkerFaceColor', 'y'); % Sole
axis equal;
xlabel('x [m]');
ylabel('y [m]');
title('Orbita della Terra attorno al Sole');
grid on;
