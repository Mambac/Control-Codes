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
ae = optParams(1:6); aed = optParams(7:12);
kmax = optParams(13:18); GE = optParams(19:24);
save('FuzzyGA.mat', 'ae', 'aed', 'kmax', 'GE');

%% ===== SIMULATION ===== %
% ===== FUZZY ===== %
fisK = FuzzyPID(ae, aed, kmax);
plotOptions = evalfisOptions(...
        'OutOfRangeInputValueMessage', 'none',...
        'NoRuleFiredMessage', 'none',...
        'EmptyOutputFuzzySetMessage', 'none');

% ===== DEFINITION OF PATH TRACKING ===== %
[qd, qdot_d, qddot_d] = PathDefined(t, Nq);

for z=1:2
    % ===== INITIAL CONDITIONS ===== %
    q = ones(Nq, length(t)); q(:,1) = [0 -pi pi pi 0 0];
    qdot = zeros(Nq, length(t));
    inte = zeros(Nq,1);
    
    taulog = zeros(Nq,length(t));
    Dlog = zeros(Nq, length(t));

    if z==1
        % ===== PID PARAMS ===== %
        load('PIDNoDist.mat', 'Kp', 'Ki', 'Kd')
    else
        % ===== FUZZY PARAMS CALCULATION ===== %
        Gu = diag(Kp)./GE;
        GCE = diag(Kd)./Gu;
        GIE = diag(Ki)./Gu;
    end

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

        if z==1
            v = qddot_d(:,i) + Kd*edot + Kp*e + Ki*inte;
        else
            u = zeros(Nq,1);
            for j = 1:Nq
                eNorm = GE(j) * e(j);
                edotNorm = GCE(j) * edot(j);
                uFuzzy = evalfis(fisK{j}, [eNorm edotNorm], plotOptions);
                u(j) = Gu(j) * ( uFuzzy + GIE(j)*inte(j) );
            end
            v = qddot_d(:,i)+u;
        end

        % === DYNAMIC COMPENSATOR === %
        tau = Mevl*v + (Cevl + Gevl);

        % === NEW STATES === %
        qddot = Mevl \ (tau - (Cevl+Gevl+F+D));
        qdot(:,i+1) = thdot + qddot*dt;
        q(:,i+1) = th + qdot(:,i+1)*dt;

        taulog(:,i) = tau; Dlog(:,i) = D;
    end

    if z==1
        qPID = q; tauPID = taulog;
    else
        qFuzzy = q; tauFuzzy = taulog;
    end
end

%% ===== ERROR ANALYSIS (after simulation) ===== %
eQ = qPID - qFuzzy;
absE = abs(eQ);

RMS  = sqrt(mean(eQ.^2, 2));
RMSMean  = mean(RMS);

RMSGlobal = sqrt(mean(eQ(:).^2));

fprintf('\n===== Error qPID vs qFuzzy =====\n');
for j = 1:Nq
    fprintf('J%d: RMS=%.6f rad', j, RMS(j));
end

fprintf('--- Joint Mean ---\n');
fprintf('RMS_mean=%.6f rad', RMSMean);
fprintf('RMS_global=%.6f rad\n\n', RMSGlobal);

%% ====== GRAPHS ====== %
titles = {'\theta_1', '\theta_2', '\theta_3', '\theta_4', '\theta_5', '\theta_6'};

figure;
sgtitle('Joints Tracking') 
for i = 1:3
    % ===== Position ===== %
    subplot(Nq/2,1,i);
    plot(t, qPID(i,:), 'b.-', t, qFuzzy(i,:), 'r--', 'LineWidth', 2);
    legend('PID','Fuzzy');
    title(['Response of ', titles{i}]); ylabel('\theta [rad]'); grid on;
    ylim([-1 1])
end
figure;
sgtitle('Joints Tracking') 
for i = 1:Nq/2
    % ===== Position ===== %
    subplot(Nq/2,1,i);
    plot(t, qPID(i+3,:), 'b.-', t, qFuzzy(i+3,:), 'r--', 'LineWidth', 2);
    legend('PID','Fuzzy');
    title(['Response of ', titles{i+3}]); ylabel('\theta [rad]'); grid on;
    ylim([-1 1])
end

figure;
sgtitle('Control signal') 
for i = 1:Nq
    % ===== Torque ===== %
    subplot(Nq,1,i);
    plot(t, tauPID(i,:), 'b-', t, tauFuzzy(i,:), 'm--', 'LineWidth', 2);
    legend('PID','Fuzzy');
    title('Torque comparison'); ylabel('\tau [Nm]'); grid on;
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

function fisK = FuzzyPID(aeVec, aedVec, kmaxVec)
    Nq = 6;
    fisK = cell(Nq, 1);
    for j = 1:Nq
        ae = aeVec(j); aed = aedVec(j); kmax = kmaxVec(j);
        fisK{j} = sugfis('Name', ['Sugeno_' num2str(j)]);

        % ===== INPUT e ===== %
        fisK{j} = addInput(fisK{j}, [-2*ae 2*ae], 'Name', 'e');
        fisK{j} = addMF(fisK{j}, 'e', 'trimf', [-3*ae -2*ae -ae], 'Name', 'NB');
        fisK{j} = addMF(fisK{j}, 'e', 'trimf', [-2*ae -ae 0], 'Name', 'N');
        fisK{j} = addMF(fisK{j}, 'e', 'trimf', [0 ae 2*ae], 'Name', 'P');
        fisK{j} = addMF(fisK{j}, 'e', 'trimf', [ae 2*ae 3*ae], 'Name', 'PB');

        % ===== INPUT edot ===== %
        fisK{j} = addInput(fisK{j}, [-2*aed 2*aed], 'Name', 'edot');
        fisK{j} = addMF(fisK{j}, 'edot', 'trimf', [-3*aed -2*aed -aed], 'Name', 'NB');
        fisK{j} = addMF(fisK{j}, 'edot', 'trimf', [-2*aed -aed 0], 'Name', 'N');
        fisK{j} = addMF(fisK{j}, 'edot', 'trimf', [0 aed 2*aed], 'Name', 'P');
        fisK{j} = addMF(fisK{j}, 'edot', 'trimf', [aed 2*aed 3*aed], 'Name', 'PB');

        % ===== OUTPUT ===== %
        fisK{j} = addOutput(fisK{j}, [-kmax kmax], 'Name', 'uFuzzy');
        fisK{j} = addMF(fisK{j}, 'uFuzzy', 'constant', -kmax, 'Name', 'NB');
        fisK{j} = addMF(fisK{j}, 'uFuzzy', 'constant', -kmax/5, 'Name', 'N');
        fisK{j} = addMF(fisK{j}, 'uFuzzy', 'constant', 0, 'Name', 'Z');
        fisK{j} = addMF(fisK{j}, 'uFuzzy', 'constant', kmax/5, 'Name', 'P');
        fisK{j} = addMF(fisK{j}, 'uFuzzy', 'constant', kmax, 'Name', 'PB');

        rules = [...
            % e = NB
            1 1 1 1 1; 1 2 1 1 1; 1 3 2 1 1; 1 4 3 1 1;
            % e = N
            2 1 1 1 1; 2 2 2 1 1; 2 3 3 1 1; 2 4 4 1 1;
            % e = P
            3 1 2 1 1; 3 2 3 1 1; 3 3 4 1 1; 3 4 5 1 1;
            % e = PB
            4 1 3 1 1; 4 2 4 1 1; 4 3 5 1 1; 4 4 5 1 1];
        
        fisK{j} = addRule(fisK{j}, rules);
    end
end