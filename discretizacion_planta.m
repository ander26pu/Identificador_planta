% identificar_planta_discreto.m
% Estima G(s)=K/(tau*s+1) de primer orden a partir de un escalón,
% usando sólo el intervalo [t0, tf], y luego discretiza con ZOH,
% Forward Euler, Backward Euler y Tustin.
% Requiere System Identification Toolbox y Control System Toolbox.

clc; clear; close all;

%% Parámetros de usuario
csvFile = 'datos_planta_primer_orden_1ms.csv';  % nombre del CSV
t0  = 0.1;    % inicio de ventana útil [s]
tf  = 7.4;    % fin de ventana útil    [s]
Ts  = 0.1;    % período de discretización [s]

%% Comprobación de toolboxes
if ~exist('iddata','file')
    error('System Identification Toolbox no está instalado.');
end
if ~exist('c2d','file')
    error('Control System Toolbox no está instalado.');
end

%% 1) Leer y preparar datos
T = readtable(csvFile,'VariableNamingRule','preserve');
if ~all(ismember({'t','u','v'}, T.Properties.VariableNames))
    error('El CSV debe tener columnas t, u, v.');
end
t = T.t;  u = T.u;  v = T.v;

% 2) Ordenar por tiempo y eliminar duplicados
[t_s,  idx1] = sort(t);
u_s = u(idx1);  v_s = v(idx1);
[t_u, idx2] = unique(t_s);
u_u = u_s(idx2);  v_u = v_s(idx2);

% 3) Crear malla uniforme con paso dt medio
dt = mean(diff(t_u));
t_all = (t_u(1):dt:t_u(end))';
u_all = interp1(t_u, u_u, t_all,'linear');
v_all = interp1(t_u, v_u, t_all,'linear');

% 4) Selección de ventana de identificación
win   = (t_all >= t0) & (t_all <= tf);
t_win = t_all(win);
u_win = u_all(win);
v_win = v_all(win);

% --- Graficar señales en la ventana ---
figure;
plot(t_all, v_all, 'Color', [0.8 0.8 0.8]); hold on;
plot(t_win, v_win, 'b', 'LineWidth', 1.5);
xlabel('Tiempo (s)'); ylabel('Salida v (V)');
title(sprintf('Ventana de Identificación: %.2f a %.2f s', t0, tf));
grid on;

%% 2) Estimación del modelo continuo de orden 1
data = iddata(v_win, u_win, dt);
sysc = tfest(data, 1, 0);

% Extraer K y tau
K   = dcgain(sysc);
p   = pole(sysc);
tau = -1/p(1);

fprintf('Modelo continuo identificado:\n');
disp(sysc);
fprintf('K = %.3f, tau = %.3f s\n', K, tau);

%% Validación continua
figure; compare(data, sysc); title('Continuo: datos vs. modelo');
figure; bode(sysc);      title('Continuo: Bode');
figure; step(sysc);      title('Continuo: Escalón');

%% 3) Convertir a espacio de estados
[Ac, Bc, Cc, Dc] = ssdata(sysc);

%% 4) Generación de sistemas discretos
sys_zoh = c2d(sysc, Ts, 'zoh');
Ad_fwd = eye(size(Ac)) + Ts*Ac;      Bd_fwd = Ts * Bc;
sys_fwd = ss(Ad_fwd, Bd_fwd, Cc, Dc, Ts);
M       = eye(size(Ac)) - Ts*Ac;     Ad_bwd = inv(M);
Bd_bwd  = Ad_bwd * Ts * Bc;
sys_bwd = ss(Ad_bwd, Bd_bwd, Cc, Dc, Ts);
sys_tus = c2d(sysc, Ts, 'tustin');

%% 5) Preparar señal de entrada y vector de tiempo discretos
t_d = (t_win(1):Ts:t_win(end))';
u_d = interp1(t_win, u_win, t_d, 'previous');

%% 6) Simulación de respuestas
y_c = lsim(sysc,    u_win, t_win);
y_z = lsim(sys_zoh, u_d,   t_d);
y_f = lsim(sys_fwd, u_d,   t_d);
y_b = lsim(sys_bwd, u_d,   t_d);
y_t = lsim(sys_tus, u_d,   t_d);

%% 7) Gráfico comparativo en un solo plot
figure; hold on; grid on;
plot(t_win, y_c, 'b-',  'LineWidth',1.5);
stairs(t_d, y_z, 'r-.','LineWidth',1);
stairs(t_d, y_f, 'm--','LineWidth',1);
stairs(t_d, y_b, 'g:', 'LineWidth',1);
stairs(t_d, y_t, 'c-', 'LineWidth',1);
xlabel('Tiempo (s)');
ylabel('Salida y');
title(sprintf('Comparación: continua vs discretas (Ts=%.3f s)', Ts));
legend({'Continua','ZOH','Euler Forward','Euler Backward','Tustin'}, ...
       'Location','best');

%% 8) Curvas 3D: y = f(t, Ts) para cada método, con color por método y brillo por Ts
Ts_vals = linspace(0.01, 1, 40);  % barrido de Ts
t_plot  = t_win;                  % usar misma malla para comparar
N_Ts    = numel(Ts_vals);
N_t     = numel(t_plot);

% Prealocar
Y_zoh = nan(N_Ts, N_t);
Y_fwd = nan(N_Ts, N_t);
Y_bwd = nan(N_Ts, N_t);
Y_tus = nan(N_Ts, N_t);

for k = 1:N_Ts
    Tsi = Ts_vals(k);
    % Discretizar
    sys_z = c2d(sysc, Tsi, 'zoh');
    AdF   = eye(size(Ac)) + Tsi*Ac;   BdF = Tsi*Bc;
    sys_f = ss(AdF, BdF, Cc, Dc, Tsi);
    M     = eye(size(Ac)) - Tsi*Ac;   AdB = inv(M); BdB = AdB*Tsi*Bc;
    sys_b = ss(AdB, BdB, Cc, Dc, Tsi);
    sys_t = c2d(sysc, Tsi, 'tustin');

    % Simulación
    td = (t_plot(1):Tsi:t_plot(end))';
    ud = interp1(t_win, u_win, td, 'previous');

    y_z = lsim(sys_z, ud, td);   Y_zoh(k,:) = interp1(td, y_z, t_plot);
    y_f = lsim(sys_f, ud, td);   Y_fwd(k,:) = interp1(td, y_f, t_plot);
    y_b = lsim(sys_b, ud, td);   Y_bwd(k,:) = interp1(td, y_b, t_plot);
    y_t = lsim(sys_t, ud, td);   Y_tus(k,:) = interp1(td, y_t, t_plot);
end

%% Colores base por método
col_z = [1 0 0];  % rojo      - ZOH
col_f = [0 0 1];  % azul      - Forward Euler
col_b = [0 0.5 0];% verde     - Backward Euler
col_t = [0 0 0];  % negro     - Tustin

%% Crear figura con líneas 3D
figure; hold on; grid on;
for k = 1:N_Ts
    alpha = 1 - (Ts_vals(k) - min(Ts_vals)) / (max(Ts_vals) - min(Ts_vals));  % 1 → oscuro, 0 → claro
    % Trazar cada curva por método
    plot3(t_plot, Ts_vals(k)*ones(size(t_plot)), Y_zoh(k,:), '-', 'Color', alpha*col_z + (1-alpha), 'LineWidth', 1);
    plot3(t_plot, Ts_vals(k)*ones(size(t_plot)), Y_fwd(k,:), '-', 'Color', alpha*col_f + (1-alpha), 'LineWidth', 1);
    plot3(t_plot, Ts_vals(k)*ones(size(t_plot)), Y_bwd(k,:), '-', 'Color', alpha*col_b + (1-alpha), 'LineWidth', 1);
    plot3(t_plot, Ts_vals(k)*ones(size(t_plot)), Y_tus(k,:), '-', 'Color', alpha*col_t + (1-alpha), 'LineWidth', 1);
end

xlabel('Tiempo t (s)');
ylabel('Periodo de muestreo Ts (s)');
zlabel('Salida y(t)');
title('Curvas 3D: y(t, Ts) por método de discretización');
view(45, 30);
legend({'ZOH','Euler Forward','Euler Backward','Tustin'}, 'Location','northeastoutside');
