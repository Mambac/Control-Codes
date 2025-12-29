clc; clear; close all;

% ===== PREPARING DYNAMICS FOR SIMULATION ====== %
fileDynamics = 'Data.mat';
if ~isfile(fileDynamics)
    disp('Calculating Dynamic of the system')
    Dynamics6();
end

load(fileDynamics);
disp('Dynamic model is ready to work!!! :D');

%% ====== PARAMETER DEFINITION ====== %
dt = 0.001;
ts = 10;
t = 0:dt:ts;
Nq = 6;
b = [0.30; 0.20; 0.12; 0.05; 0.05; 0.02]; % Viscous constants (Nm s/rad)
c = [1.50; 1.00; 0.60; 0.20; 0.02; 0.02]; % Coulomb constants (Nm)

% ===== MATRICES ===== %
Mq_func = matlabFunction(M, 'Vars', {sym('q',[Nq 1],'real')});
Cq_func = matlabFunction(C,'Vars',{sym('q',[Nq 1],'real'), ...
                         sym('dq',[Nq 1],'real')},'Optimize',false);
Gq_func = matlabFunction(G, 'Vars', {sym('q',[Nq 1],'real')});

%% ===== OPTIMIZED BEHAVIOR ===== %

%{
    These params supposed be optimized using any optimizer of your choice.
    ex: GA, MOGA, PSO, GWO, etc...
%}
Kp = diag(optParams(1:6));
Ki = diag(optParams(7:12));
Kd = diag(optParams(13:18));
save('PIDNoDist.mat', 'Kp', 'Ki', 'Kd')

% ===== INITIAL CONDITIONS ===== %
q = ones(Nq, length(t)); q(:,1) = [0 -pi pi pi 0 0];
qdot = zeros(Nq, length(t));
inte = zeros(Nq,1);

% ===== DEFINITION OF PATH TRACKING ===== %
[qd, qdot_d, qddot_d] = PathDefined(t, Nq);

Slog = zeros(Nq,length(t)); SsatLog = zeros(Nq,length(t)); taulog = zeros(Nq,length(t));
Dlog = zeros(Nq, length(t));


% ===== MAIN LOOP ===== %
for i=1:length(t)-1
    th = q(:,i); thdot = qdot(:,i);

    % === DYNAMIC DEFINITION === %
    Mevl = Mq_func(th);
    Cmat = Cq_func(th,thdot);
    if any(isnan(Cmat(:)))
        disp('NaN in Cmat at t = ' + string(t(i)));
        break;
    end
    Cevl = Cmat*thdot;
    Gevl = Gq_func(th);

    % === FRICTION AND DISTURBANCES === %
    F = b.*thdot + c.*sign(thdot);
    D = [0; 0; 0; 0; 0; 0];

    % === CONTROL === %
    e = qd(:,i) - th;
    edot = qdot_d(:,i) - thdot;
    inte = inte + e*dt;

    v = qddot_d(:,i) + Kd*edot + Kp*e + Ki*inte;
    tau = Mevl*v + (Cevl + Gevl);

    % === NEW STATES === %
    qddot = Mevl \ (tau - (Cevl+Gevl+F+D));
    qdot(:,i+1) = thdot + qddot*dt;
    q(:,i+1) = th + qdot(:,i+1)*dt;

    taulog(:,i) = tau; Dlog(:,i) = D; 
end

%% ====== GRAPHS ====== %
titles = {'\theta_1', '\theta_2', '\theta_3', '\theta_4', '\theta_5', '\theta_6'};

for i = 1:Nq
    figure;
    sgtitle(['Joint ', num2str(i)])
    
    % ===== Position ===== %
    subplot(3,1,1);
    plot(t, q(i,:), 'r--', 'LineWidth', 2); hold on;
    plot(t, qd(i,:), 'b:', 'LineWidth', 2);
    legend('Output', 'Reference', 'Location','best');
    title(['Answer of ', titles{i}]);
    ylabel('\theta [rad]');
    ylim([-1 1])
    grid on;
    
    % ===== Disturbances ===== %
    subplot(3,1,2);
    plot(t, Dlog(i,:), 'g--', 'LineWidth', 1.5);
    title('Disturbance D');
    ylabel('Disturbance torque [Nm]');
    grid on;
    
    % ===== Torque ===== %
    subplot(3,1,3);
    plot(t, taulog(i,:), 'm', 'LineWidth', 1.5);
    title('Control torque \tau');
    xlabel('Time [s]');
    ylabel('\tau [Nm]');
    grid on;
end

%% ===== FUNCTIONS ====== %
function [qd, qdot_d, qddot_d] = PathDefined(t, Nq, varargin)

    p = inputParser;
    addParameter(p, 'seed', randi(1000), @isscalar);
    parse(p, varargin{:});
    rng(p.Results.seed);

    dt = t(2) - t(1);
    N = length(t);
    qd = zeros(Nq, N);

    for j = 1:Nq
        n_puntos = 12 + randi(5);
        tKey = sort(rand(1, n_puntos) * 9.5) + 0.25;
        tKey = [0, tKey, 10];

        
        amp = 0.18 + 0.1*rand;
        qKey = amp * (2*rand(1, length(tKey)) - 1);
        qKey(1) = -0.2; qKey(end) = -0.2;
        qKey = qKey - 0.2;

        
        qd(j,:) = interp1(tKey, qKey, t, 'pchip')';

        
        if j > 1
            shift = (j-1)*0.6 + 0.3*rand;
            tShift = t - shift;
            tShift(tShift < 0) = 0;
            tShift(tShift > 10) = 10;
            q_temp = interp1(tKey, qKey, tShift, 'pchip', 'extrap')';
            qd(j,:) = q_temp;
        end
    end

    
    qdot_d = [zeros(Nq,1), diff(qd,1,2)/dt];
    qddot_d = [zeros(Nq,2), diff(qdot_d,1,2)/dt];

    for j = 1:Nq
        qddot_d(j,:) = smoothdata(qddot_d(j,:), 'movmean', 11);
    end
end
