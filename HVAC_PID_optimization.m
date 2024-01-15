clear all; close all; clc;
set(0,'defaultTextInterpreter','latex')
warning('on', 'globaloptim:particleswarm:initialSwarmNotInBounds');
set(0,'defaultTextInterpreter','latex')

global T_o_data T_z_data Q_P M w dw

HVAC_par %Load model
Q_P = Q_P_data;

w = 295*ones(10000,1);
w(2000:5000) = 296;
w(8000:10000) = 294;
dw = 2;
wlb = w - 2;
wub = w + 2;

if 0

    %Fixed parameters
    par_0 = [M.u_1, M.u_2, M.u_3];

    %Constraints
    u_1_min = -1e0;
    u_1_max = -1e-6;
    u_2_min = -1e-1;
    u_2_max = -1e-6;
    u_3_min = 2e-3;
    u_3_max = 1e1;

    lb = [u_1_min, u_2_min, u_3_min];
    ub = [u_1_max, u_2_max, u_3_max];

    options = optimoptions('particleswarm','maxIter', 20, 'Display','iter','SwarmSize',100,'PlotFcn', @custom_pswplotbestf);
    options.InitialSwarmMatrix = par_0;
    [par, fval, exitflag, e] = particleswarm(@func, length(par_0), lb, ...
        ub, options);

    save('HVAC_PID','par')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HVAC_par

load HVAC_PID
M.u_1 = par(1);
M.u_2 = par(2);
M.u_3 = par(3);

k_0 = max([M.d_I,M.d_D])+2;
u = 0*ones(k_0,1);
Q_I = 0*ones(k_0,1);
Q_D = 0*ones(k_0,1);
T_z = T_z_data(1:k_0);
sum_T_z = 0*ones(k_0,1);
T_z_0 = T_z_data;
u_sum = 0*ones(k_0,1);

for k = k_0:1:length(w)

    sum_T_z(k) = sum_T_z(k-1) + (w(k)-T_z(k-1));
    u(k) = M.u_1*(w(k)-T_z(k-1))+  M.u_2*sum_T_z(k);
    [T_z(k), Q_I(k), Q_D(k)] = HVAC(T_z(k-1), T_z(k-M.d_I-1), T_z(k-M.d_D-1),...
        u(k-1-M.d_I), T_o_data(k-M.d_I-1), T_o_data(k-M.d_D-1), Q_I(k-1), ...
        Q_D(k-1) , Q_P(k-1),  M);
    u_sum(k) = u_sum(k-1) + u(k)^2;

end


h1 = figure(1);
subplot(4,1,1:2); hold off;
p1 = plot(T_z,'r'); hold on;
p2 = plot(wub,'k');
plot(wlb,'k')
p3 = plot(w,'k--');
ax = gca;
ax.XAxis.Exponent = 4;
%     hcaX=ax.XRuler;
%     hcaX.SecondaryLabel.Position(2) = hcaX.SecondaryLabel.Position(2)*1.04;
%     hcaX.SecondaryLabel.Position(1) = hcaX.SecondaryLabel.Position(1)*1.04;
xlim('tight')
ylabel('Temperature')
legend([p1,p2,p3],'System output', 'Reference interval','Reference value','Location','best')

subplot(4,1,3); hold off; hold off;
h2 = plot(u,'r');
ax = gca;
ax.XAxis.Exponent = 4;
%     hcaX=ax.XRuler;
%     hcaX.SecondaryLabel.Position(2) = hcaX.SecondaryLabel.Position(2)*1.024;
%     hcaX.SecondaryLabel.Position(1) = hcaX.SecondaryLabel.Position(1)*1.04;
ylabel({'Control','signal'})
xlim('tight')

subplot(4,1,4); hold off;
plot(u_sum,"color",'r');
text(length(u_sum),u_sum(end),['$\leftarrow$',num2str(ceil(u_sum(end)))],'FontSize',8)
ax = gca;
ax.XAxis.Exponent = 4;
ylabel('Cost')
legend('Sum of squared inputs','Location','best')
xlabel('Time step')
xlim('tight')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% OPTIMIZATION FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function err = func(par)

% global K_D_min K_D_max tau_D_min tau_D_max  K_I_min ...
%     K_I_max tau_I_min tau_I_max a_I_min a_I_max b_I_min b_I_max
global T_o_data T_z_data Q_P M w dw

wlb = w - 2;
wub = w + 2;

M.u_1 = par(1);
M.u_2 = par(2);
M.u_3 = par(3);


k_0 = max([M.d_I,M.d_D])+2;
u = 0*ones(k_0,1);
Q_I = 0*ones(k_0,1);
Q_D = 0*ones(k_0,1);
T_z = T_z_data(1:k_0);
sum_T_z = 0*ones(k_0,1);
u_sum = 0;
for k = k_0:1:length(w)

    sum_T_z(k) = sum_T_z(k-1) + (w(k)-T_z(k-1));
    u(k) = M.u_1*(w(k)-T_z(k-1))+  M.u_2*sum_T_z(k);
    u_sum = u_sum + u(k)^2;

    [T_z(k), Q_I(k), Q_D(k)] = HVAC(T_z(k-1), T_z(k-M.d_I-1), T_z(k-M.d_D-1),...
        u(k-1-M.d_I), T_o_data(k-M.d_I-1), T_o_data(k-M.d_D-1), Q_I(k-1),...
        Q_D(k-1) , Q_P(k-1),  M);

end

err = u_sum(end);

cons = any(T_z <= wlb) || any(T_z >= wub);
if cons
    err = err + 1e10;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% DISPLAY FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function stop = custom_pswplotbestf(optimValues, state)

stop = false; % This function should not stop the algorithm

switch state
    case 'init'
    case 'iter'
        global T_o_data T_z_data Q_P M w dw

        wlb = w - 2;
        wub = w + 2;

        M.u_1 = optimValues.bestx(1);
        M.u_2 = optimValues.bestx(2);
        M.u_3 = optimValues.bestx(3);

        k_0 = max([M.d_I,M.d_D])+2;
        u = 0*ones(k_0,1);
        Q_I = 0*ones(k_0,1);
        Q_D = 0*ones(k_0,1);
        T_z = T_z_data(1:k_0);
        sum_T_z = 0*ones(k_0,1);
        u_sum = 0*ones(k_0,1);

        for k = k_0:1:length(w)

            sum_T_z(k) = sum_T_z(k-1) + (w(k)-T_z(k-1));
            u(k) = M.u_1*(w(k)-T_z(k-1))+  M.u_2*sum_T_z(k);

            [T_z(k), Q_I(k), Q_D(k)] = HVAC(T_z(k-1), T_z(k-M.d_I-1), T_z(k-M.d_D-1),...
                u(k-1-M.d_I), T_o_data(k-M.d_I-1), T_o_data(k-M.d_D-1), Q_I(k-1),...
                Q_D(k-1) , Q_P(k-1), M);
            u_sum(k) = u_sum(k-1) + u(k)^2;

        end

        h1 = figure(1);
        subplot(4,1,1:2); hold off;
        p1 = plot(T_z,'r'); hold on;
        p2 = plot(wub,'k');
        plot(wlb,'k')
        p3 = plot(w,'k--');
        ax = gca;
        ax.XAxis.Exponent = 4;
        %     hcaX=ax.XRuler;
        %     hcaX.SecondaryLabel.Position(2) = hcaX.SecondaryLabel.Position(2)*1.04;
        %     hcaX.SecondaryLabel.Position(1) = hcaX.SecondaryLabel.Position(1)*1.04;
        xlim('tight')
        ylabel('Temperature')
        legend([p1,p2,p3],'System output', 'Reference interval','Reference value','Location','best')

        subplot(4,1,3); hold off; hold off;
        h2 = plot(u,'r');
        ax = gca;
        ax.XAxis.Exponent = 4;
        %     hcaX=ax.XRuler;
        %     hcaX.SecondaryLabel.Position(2) = hcaX.SecondaryLabel.Position(2)*1.024;
        %     hcaX.SecondaryLabel.Position(1) = hcaX.SecondaryLabel.Position(1)*1.04;
        ylabel({'Control','signal'})
        xlim('tight')

        subplot(4,1,4); hold off;
        plot(u_sum,"color",'r');
        text(length(u_sum),u_sum(end),['$\leftarrow$',num2str(ceil(u_sum(end)))],'FontSize',8)
        ax = gca;
        ax.XAxis.Exponent = 4;
        ylabel('Cost')
        legend('Sum of squared inputs','Location','best')
        xlabel('Time step')
        xlim('tight')

        name = '..\HVAC_Images\HVAC_PID_optimal.pdf';
        exportgraphics(gcf,name,'BackgroundColor','none');
        pause(0.0001)

    case 'done'
        % Finalize your plot here
        hold off;
end

end
