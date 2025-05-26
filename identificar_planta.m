% identificar_planta.m (v3)
% Estima un sistema de primer orden G(s)=K/(tau*s+1)
% a partir de un experimento de escalón único,
% usando solo el intervalo [t0, tf] para la identificación.
% Requiere System Identification Toolbox.

clc; clear; close all;

% --- Comprueba toolbox ---
if ~exist('iddata','file')
    error('System Identification Toolbox no está instalado.');
end

% --- Parámetros ---
csvFile = 'datos_planta_segundo_orden_1ms.csv';
% Tiempo inicial y final (s) para identificación:
t0 = 0.0;    % segundos: inicio de la ventana útil
tf = 17.0;    % segundos: fin de la ventana útil

% --- Leer datos ---
opts = detectImportOptions(csvFile, 'VariableNamingRule','preserve');
T = readtable(csvFile, opts);

% Verifica columnas
req = {'t','u','v'};
if ~all(ismember(req, T.Properties.VariableNames))
    error('El CSV debe tener encabezados: t, u, v');
end

t = T.t;     % tiempo en segundos
u = T.u;     % entrada (escalón)
v = T.v;     % salida medida

% --- Ordenar y eliminar duplicados de tiempo ---
[t_sorted, sortIdx] = sort(t);
u_sorted = u(sortIdx);
v_sorted = v(sortIdx);
[t_unique, uniqueIdx] = unique(t_sorted);
u_unique = u_sorted(uniqueIdx);
v_unique = v_sorted(uniqueIdx);

% --- Definir muestreo uniforme ---
dt = mean(diff(t_unique));
t_uniform = (t_unique(1):dt:t_unique(end))';

% --- Interpolación sobre tiempos únicos ---
u_uniform = interp1(t_unique, u_unique, t_uniform, 'linear');
v_uniform = interp1(t_unique, v_unique, t_uniform, 'linear');

% --- Selección de ventana de identificación ---
idx_win = (t_uniform >= t0) & (t_uniform <= tf);
t_win = t_uniform(idx_win);
u_win = u_uniform(idx_win);
v_win = v_uniform(idx_win);

% --- Crear iddata solo con la ventana seleccionada ---
data = iddata(v_win, u_win, dt);

% --- Graficar señales en la ventana ---
figure;
plot(t_uniform, v_uniform, 'Color', [0.8 0.8 0.8]); hold on;
plot(t_win, v_win, 'b', 'LineWidth', 1.5);
xlabel('Tiempo (s)'); ylabel('Salida v (V)');
title(sprintf('Ventana de Identificación: %.2f a %.2f s', t0, tf));
grid on;

% --- Estimación de orden n ---
sys = tfest(data, 2, 0);

% --- Resultados ---
disp('Modelo identificado:');
disp(sys);
K   = dcgain(sys);
tau = -1/pole(sys);
fprintf('Ganancia DC K  = %.3f\n', K);
fprintf('Constante tau  = %.3f s\n', tau);

% --- Validación ---
figure; compare(data, sys);
title('Datos vs. Modelo Identificado');

figure; bode(sys);
title('Bode del Modelo');

figure; step(sys);
title('Respuesta al Escalón Modelo');
