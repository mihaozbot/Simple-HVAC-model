clear all; close all; clc;
set(0,'defaultTextInterpreter','latex')
warning('on', 'globaloptim:particleswarm:initialSwarmNotInBounds');
set(0,'defaultTextInterpreter','latex')

global ts c_z % d_D d_I
%global K_D_min K_D_max tau_D_min tau_D_max K_I_min ...
%        K_I_max tau_I_min tau_I_max a_I_min a_I_max b_I_min b_I_max ...
%        u_1_min u_2_min u_3_min d_D_min d_I_min
global T_o_data T_z_data u_data Q_P_0 T_0

load('data_HVAC') %T_o_data T_z_data u_data

Q_P_0 = (0.5+0.5*sin((90:1:length(u_data(:,1))+90-1)*(2*pi)/144)); %+ ... 0.5*sin((110:1:length(u_data(:,1))+110-1)*(2*pi)/10080));
% Q_P_weekend = [3*1440:10080:length(Q_P);5*1440:10080:length(Q_P)]';
% for iw = 1:1:length(Q_P_weekend(:,1))
%     Q_P(Q_P_weekend(iw,1):Q_P_weekend(iw,2)) = 0;
% end
%
Q_P_week = [15:144:(length(Q_P_0)-144);(144+15):144:length(Q_P_0)]';

%Q_P = u_3*((T_z_data-u_0')'+ abs(min(T_z_data-u_0')'));
% figure(11); hold off;
% plot(T_z_data(Q_P_week(1,1):Q_P_week(1,2)));  hold on;
% plot(min(T_z_data(Q_P_week(1,1):Q_P_week(1,2))) + Q_P(Q_P_week(1,1):Q_P_week(1,2))/u_3);
T_0 = 295*ones(1,length(u_data(:,1)));

for iw = 1:1:length(Q_P_week(:,1))
    % if (sum((T_z_data(Q_P_week(iw,1):Q_P_week(iw,2)) - min(T_z_data(Q_P_week(iw,1):Q_P_week(iw,2))) - ...
    %         Q_P_0(Q_P_week(iw,1):Q_P_week(iw,2))').^2)> ...
    %         sum((T_z_data(Q_P_week(iw,1):Q_P_week(iw,2)) - min(T_z_data(Q_P_week(iw,1):Q_P_week(iw,2)))).^2))
    %     Q_P_0(Q_P_week(iw,1):Q_P_week(iw,2)) = 0;
    % else
    T_o_0_test(Q_P_week(iw,1):Q_P_week(iw,2)) = (max(abs(T_o_data(Q_P_week(iw,1):Q_P_week(iw,2))-T_z_data(Q_P_week(iw,1):Q_P_week(iw,2))))...
        - max(abs(T_o_data(Q_P_week(iw,1))-T_z_data(Q_P_week(iw,1))),abs(T_o_data(Q_P_week(iw,2))-T_z_data(Q_P_week(iw,2)))));
    Q_P_0_test(Q_P_week(iw,1):Q_P_week(iw,2)) = (max(T_z_data(Q_P_week(iw,1):Q_P_week(iw,2))) - max(T_z_data(Q_P_week(iw,1)),T_z_data(Q_P_week(iw,2))));
    Q_P_0(Q_P_week(iw,1):Q_P_week(iw,2)) = Q_P_0(Q_P_week(iw,1):Q_P_week(iw,2))*...
        (min(3,max(T_z_data(Q_P_week(iw,1):Q_P_week(iw,2))) - max(T_z_data(Q_P_week(iw,1)),T_z_data(Q_P_week(iw,2)))));
    % end
    T_0(Q_P_week(iw,1):Q_P_week(iw,2)) = max(294, min(300, mean(T_z_data(Q_P_week(iw,1):Q_P_week(iw,2))) ) );
end
T_0(18591:end) = T_0(18591:end);

figure(11); hold off;
plot(T_z_data(1:end)); hold on;
plot(295 + Q_P_0(1:end));

%Fixed parameters
ts = 60*10;
%k_z = 20;
c_z = 72e4;


if 0
    %Optimized parameters
    if 0
        K_D = 0.7;
        tau_D = 600.0588;
        K_I = 1.8166;
        tau_I = 600.2028;
        a_I =  -6.0609 ;
        b_I =  -22.3376;
        u_1 = -0.26 ;
        u_2 =  -0.0001;
        u_3 =  0.9454;
        d_D = 1;
        d_I = 2;
    else
        load('HVAC_par.mat')
        K_D = par(1);
        tau_D = par(2);
        K_I = par(3);
        tau_I = par(4);
        a_I = par(5);
        b_I = par(6);
        u_1 = par(7);
        u_2 = par(8);
        u_3 = par(9);
        d_D = round(par(10));
        d_I = round(par(11));

    end

    par_0 = [K_D, tau_D, K_I, tau_I, a_I, b_I, u_1, u_2, u_3, d_D, d_I];
    %-0.0147    1.5324    1.8761    1.3569    0.4108    1.6934    1.1552   -1.1223    1.4481    0.7957    3.4114

    %Constraints
    K_D_min = 6e0;
    K_D_max = 6e0;
    tau_D_min = 10*ts;
    tau_D_max = 30*ts;
    K_I_min = 3e-2;
    K_I_max = 4e0;
    tau_I_min = 10*ts;
    tau_I_max = 30*ts;
    a_I_min = -3e0;
    a_I_max = 3e0;
    b_I_min = -5e2;
    b_I_max = 5e2;
    u_1_min = -1e0;
    u_1_max = -1e-6;
    u_2_min = -1e-1;
    u_2_max = -1e-6;
    u_3_min = 2e-3;
    u_3_max = 1e1;
    d_D_min = 1e0;
    d_D_max = 1e1;
    d_I_min = 2e0;
    d_I_max = 2e1;

    lb = [K_D_min, tau_D_min, K_I_min, tau_I_min, a_I_min, b_I_min, u_1_min,...
        u_2_min, u_3_min, d_D_min, d_I_min];
    ub = [K_D_max, tau_D_max, K_I_max, tau_I_max, a_I_max, b_I_max, u_1_max,...
        u_2_max, u_3_max, d_D_max, d_I_max];

    %     options = optimset('maxIter', 100,'Display','iter');
    %     [par, fval, exitflag, e] = fminsearch(@func, par_0, options);
    options = optimoptions('particleswarm','maxIter', 200, 'Display','iter','SwarmSize',100,'PlotFcn', @custom_pswplotbestf);%);

    options.InitialSwarmMatrix = par;
    [par, fval, exitflag, e] = particleswarm(@func, length(par_0), lb, ...
        ub, options);

    save('HVAC_par','par')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HVAC_par

k_0 = max([M.d_I,M.d_D])+2;
u = 0*ones(k_0,1);
Q_I = 0*ones(k_0,1);
Q_D = 0*ones(k_0,1);
T_z = T_z_data(1:k_0);
k_min = 1000;
sum_T_z = 0*ones(k_0,1);
%u = (u_1.*u_data(:,1) + u_2.*u_data(:,2)).*u_data(:,3);
%u = u_data;
% T_0 = 295*ones(1,length(u_data(:,1)));
% if 1
%     T_0(3566:4802) = 294.75;
%     % %u_0(48023:50330) = 295.5;
%     T_0(5839:6143) = 296.5;
%     T_0(6138:15425) = 296;
%     %u_0(122125:137650) = 296.6;
%     %T_0(16535:17059) = 295;
%     %u_0(15425:18125) = 295.5;
%     T_0(15425:end) = 295;
%
% else
%     T_0(5805:15431) = 296;
% end
T_z_0 = T_z_data;
Q_P = Q_P_data;

for k = k_0:1:(length(u_data(:,1))-M.d_I)


    %sum_T_z(k) = sum_T_z(k-1) + (T_z_data(k-1+d_I)-T_z(k-1));
    %u(k) = a_I*(u_1*(T_z_data(k-1+d_I+round(tau_I/ts))-T_z(k-1)) + u_2*sum_T_z(k-1));
    %max_u = 0.1;
    %u(k-1) = u(k-2) + max(min(u(k-1)-u(k-2),max_u),-max_u);
    %sum_T_z(k) = sum_T_z(k-1) + (u_0(k)-T_z(k-1));
    %u(k-1) = u_3*(u_0(k)-T_z(k-1)) + u_2*sum_T_z(k) + u_1*(T_z_data(k-1+d_I)-T_z(k-1));
    %     if abs(T_0(k)-T_z(k-1))>0.5
    sum_T_z(k) = sum_T_z(k-1) + (T_0(k)-T_z(k-1));
    u(k) = M.u_1*(T_0(k)-T_z(k-1))+  M.u_2*sum_T_z(k);
    %     else
    %         sum_T_z(k) = 0;
    %         u(k) = 0;
    %     end
    [T_z(k), Q_I(k), Q_D(k)] = HVAC(T_z(k-1), T_z(k-M.d_I-1), T_z(k-M.d_D-1),...
        u(k-1-M.d_I), T_o_data(k-M.d_I-1), T_o_data(k-M.d_D-1),...
        Q_I(k-1), Q_D(k-1) , Q_P(k-1),  M);

end

if 1

    subplot(4,1,1:2); hold off;
    plot(T_z(1:end),'r');hold on;
    plot(T_z_data(1:length(T_z)),'b')
    xlim('tight')
    legend('Room temperature model','Room temperature measurements','Location','best')
    ylabel({'Temperature [°K]'});
    hca=gca;
    hcaX=hca.XRuler;
    hcaX.SecondaryLabel.Position(2) = hcaX.SecondaryLabel.Position(2)*1.004;
    hcaX.SecondaryLabel.Position(1) = hcaX.SecondaryLabel.Position(1)*1.04;

    subplot(4,1,3); hold off;
    plot(Q_I(1:end),'r'); hold on;
    plot(Q_D(1:end),'b')
    plot(Q_P(1:end),'k')
    ylabel({'Thermal','energy'})
    xlim('tight')
    legend('Energy controler','Energy outside','Energy people','Location','best')
    subplot(4,1,4); hold off;
    plot(u(1:end),'b'); hold on;
    ylabel({'Control','signal'}); xlabel('Time step')
    xlim('tight')
    name = '..\HVAC_Images\HVAC_model_with_PID.pdf';
    exportgraphics(gcf,name,'BackgroundColor','none');
    pause(0.0001)

end

segments = [2200,2800; 4600,5200; 9500,10000; 21000, 22000];
len = length(T_z);
figure(31);
for i = 1:4
    k_seg = segments(i,1):1:segments(i,2);
    subplot(2, 2, i); hold off;
    plot(k_seg, T_z(k_seg),'r'); hold on;
    plot(k_seg, T_z_data(k_seg),'b');
    xlabel('Time step')
    xlim('tight')

    ax = gca;
    ax.XAxis.Exponent = 4;
    if i ==2
        legend('Room temperature model','Room temperature measurements','Location','best')
    end
    name = '..\HVAC_Images\HVAC_model_zoom_in.pdf';
    exportgraphics(gcf,name,'BackgroundColor','none');
end

figure(32)
subplot(2,1,1);
plot(u,atan(u),'.r'); hold on; xlim("tight"); ylim("tight"); axis manual;
plot([-2,2],[-2,2],'b');
legend('Nonlinearity arctan($u$)','Linear function','Location','best','interpreter','latex')
name = '..\HVAC_Images\HVAC_model_input.pdf';
exportgraphics(gcf,name,'BackgroundColor','none');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k_0 = max([M.d_I,M.d_D])+2;
u = 0*ones(k_0,1);
Q_I = 0*ones(k_0,1);
Q_D = 0*ones(k_0,1);
T_z = T_z_data(1:k_0);
k_min = 1000;

for k = k_0:1:length(u_data)-M.d_I
    u(k-1) = 0;

    %max_u = 0.1;
    %u(k-1) = u(k-2) + max(min(u(k-1)-u(k-2),max_u),-max_u);
    %sum_T_z(k) = sum_T_z(k-1) + (T_z(k-1)-u_0(k));
    %u(k-1-d_I) = u_1*(T_z(k-1)-u_0(k)) + u_2*sum_T_z(k);
    [T_z(k), Q_I(k), Q_D(k)] = HVAC(T_z(k-1), T_z(k-M.d_I-1), T_z(k-M.d_D-1), u(k-1-M.d_I), T_o_data(k-M.d_I-1), T_o_data(k-M.d_D-1), Q_I(k-1), Q_D(k-1) , 0*Q_P(k-1),  M);

end

if 1
    figure(3);
    subplot(4,1,1:2); hold off;
    plot(T_o_data); hold on;
    plot(T_z(1:end),'r');
    plot(T_z_data(1:length(T_z)),'b')
    xlim('tight')
    legend('Outside temperature measurements','Room temperature model','Room temperature measurements','Location','best')
    ylabel({'Temperature'}); xlabel('Time step')
    hca=gca;
    hcaX=hca.XRuler;
    hcaX.SecondaryLabel.Position(2) = hcaX.SecondaryLabel.Position(2)*1.024;
    hcaX.SecondaryLabel.Position(1) = hcaX.SecondaryLabel.Position(1)*1.04;
    subplot(4,1,3); hold off;
    plot(Q_I(1:end),'r'); hold on;
    plot(Q_D(1:end),'b')
    plot(Q_P(1:end),'k')
    ylabel({'Thermal','energy'})
    xlim('tight')
    legend('Energy controler','Energy outside','Energy people','Location','best')
    name = '..\HVAC_Images\HVAC_model_outside_temperature.pdf';
    exportgraphics(gcf,name,'BackgroundColor','none');
    pause(0.0001)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% k_0 = max([M.d_I,M.d_D])+2;
% u = 0*ones(k_0,1);
% Q_I = 0*ones(k_0,1);
% Q_D = 0*ones(k_0,1);
% T_z = T_z_data(1:k_0);
% sum_T_z = 0*ones(k_0,1);
%
% if 1
%     u_0(35686:48023) = 294.75;
%     % %u_0(48023:50330) = 295.5;
%     u_0(58394:61431) = 296.5;
%     u_0(61385:154250) = 296;
%     u_0(73293:82851) = 296;
%     %u_0(122125:137650) = 296.6;
%     u_0(165350:170590) = 295;
%     u_0(154250:181253) = 295.5;
%     %u_0(181253:end) = 295;
% else
%     u_0(58205:154311) = 296;
% end
% T_z_0 = T_z_data;
% Q_P = M.u_3*(0.5+0.5*sin((110:1:length(u_data(:,1))+110-1)*(2*pi)/144));% ...
% %0.5*sin((110:1:length(u_data(:,1))+110-1)*(2*pi)/1008));
%
% for k = k_0:1:length(u_data)-M.d_I
%
%     sum_T_z(k) = sum_T_z(k-1) + (T_z(k-1)-u_0(k-1));
%     u(k) = 0.1*(u_0(k-1)-T_z(k-1)) + 0.001*sum_T_z(k);
%
%     [T_z(k), Q_I(k), Q_D(k)] = HVAC(T_z(k-1), T_z(k-M.d_D-1), u(k-1-M.d_I), T_o_data(k-M.d_I-1), Q_I(k-1), Q_D(k-1) , Q_P(k),  M);
%
% end
%
% if 1
%     figure(4);
%     subplot(4,1,1:2); hold off;
%     plot(T_z(1:end),'r'); hold on;
%     plot(T_z_data(1:end),'b');
%     plot(u_0,'k')
%     legend('Room temperature model','Room temperature measurements','Reference','Location','best')
%     xlim('tight')
%     subplot(4,1,3); hold off;
%     plot(Q_I(1:end),'r'); hold on;
%     plot(Q_D(1:end),'b')
%     xlim('tight')
%     legend('Energy controler','Energy outside','Location','best')
%     subplot(4,1,4); hold off;
%     plot(u(1:end),'b'); hold on;
%     legend('Control signal','Location','best')
%     xlim('tight')
%     ylabel('$u$'); xlabel('Time step')
%     pause(0.0001)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% OPTIMIZATION FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function err = func(par)

global ts k_z c_z% d_D d_I
% global K_D_min K_D_max tau_D_min tau_D_max  K_I_min ...
%     K_I_max tau_I_min tau_I_max a_I_min a_I_max b_I_min b_I_max
global T_o_data T_z_data u_data Q_P_0 T_0

K_D = par(1);
tau_D = par(2);
K_I = par(3);
tau_I = par(4);
a_I = par(5);
b_I = par(6);
u_1 = par(7);
u_2 = par(8);
u_3 = par(9);
d_D = abs(round(par(10)));
d_I = abs(round(par(11)));

k_0 = max([d_I,d_D])+2; % +0.7*length(u_data(:,1));
u = 0*ones(k_0,1);
Q_I = 0*ones(k_0,1);
Q_D = 0*ones(k_0,1);
T_z = T_z_data(1:k_0);
k_min = 1000;
sum_T_z = 0*ones(k_0,1);
K_O = 0*ones(k_0,1)>1;

%u = (u_1.*u_data(:,1) + u_2.*u_data(:,2)).*u_data(:,3);
%u = u_data;
% T_0 = 295*ones(1,length(u_data(:,1)));
% if 1
%     T_0(3566:4802) = 294.75;
%     % %u_0(48023:50330) = 295.5;
%     T_0(5839:6143) = 296.5;
%     T_0(6138:15425) = 296;
%     %u_0(122125:137650) = 296.6;
%     %u_0(16535:17059) = 295;
%     T_0(15425:end) = 295;
%
% else
%     T_0(5805:15431) = 295.5;
% end

T_z_0 = max(294, min(300, T_z_data));
u_th = 0.9;
%T_z_0(T_z_0>(u_0'+u_th)) = u_0(T_z_0>(u_0'+u_th))'+u_th; % - 0.5*randn(length(T_z_data), 1);
%T_z_0(T_z_0<(u_0'-u_th)) = u_0(T_z_0<(u_0'-u_th))'-u_th; % - 0.5*randn(length(T_z_data), 1);
Q_P = u_3*Q_P_0;

for k = k_0:1:(length(u_data(:,1))-d_I-round(tau_I/ts))

    %sum_T_z(k) = sum_T_z(k-1) + (T_z_data(k-1+d_I)-T_z(k-1));
    %u(k) = a_I*(u_1*(T_z_data(k-1+d_I+round(tau_I/ts))-T_z(k-1)) + u_2*sum_T_z(k-1)); %TLE MROE BIT U(K-1)
    %max_u = 0.1;
    %u(k-1) = u(k-2) + max(min(u(k-1)-u(k-2),max_u),-max_u);

    %     if abs(T_0(k)-T_z(k-1))>0.5
    sum_T_z(k) = sum_T_z(k-1) + (T_0(k)-T_z(k-1));
    u(k) = u_1*((T_0(k)-T_z(k-1)) + u_2*sum_T_z(k));
    %u(k-d_I) = u_1*(T_0(k)-T_z(k-1)) + u_2*(T_z_data(k-1)-T_z(k-1));
    %     else
    %         sum_T_z(k) = 0;
    %         u(k) = 0;
    %     end
    %K_O(k) = ((T_o_data(k-1) > T_z(k-1)) && (T_0(k-1)>T_z(k-1))) || ((T_o_data(k-1) < T_z(k-1)) && (T_0(k-1)<T_z(k-1)));
    Q_I(k) = Q_I(k-1) + ts*(K_I*atan(u(k-1-d_I))*(b_I + (T_o_data(k-d_I-1)-T_z(k-d_I-1)))-Q_I(k-1))/tau_I; %Energy input
    Q_D(k) = Q_D(k-1) + ts*(K_D*(T_o_data(k-d_D-1)-T_z(k-d_D-1))-Q_D(k-1))/tau_D; %Energy disturbance
    T_z(k) = T_z(k-1) + ts*(Q_I(k) + Q_D(k) + Q_P(k))/c_z;  %Zone temperature

end
% err_T_z = max([sum((T_z(k_min:k)-T_z_0(k_min:k)- 0.1*randn(k-k_min+1, 1)).^2),...
%     sum((T_z(k_min:k)-T_z_0(k_min:k)- 0.1*randn(k-k_min+1, 1)).^2),...
%     sum((T_z(k_min:k)-T_z_0(k_min:k)- 0.05*randn(k-k_min+1, 1)).^2),...
%     sum((T_z(k_min:k)-T_z_0(k_min:k)- 0.05*randn(k-k_min+1, 1)).^2)]);

T = ts;   % Sampling period
Fs = 1/ts;            % Sampling frequency
L = (length(u_data(:,1))-d_I-round(tau_I/ts));             % Length of signal
t = (0:(L-1))*T;        % Time vector
T_z_0_fft = fft(T_z_0(k_min:k));
T_z_fft = fft(T_z(k_min:k));

if 0
    figure(42); hold off;
    plot((abs(T_z_0_fft(1:end/2)))); hold on;
    plot((abs(T_z_fft(1:end/2))));
end

%  real( (T_z_0_fft(10:floor(end/2)))-(T_z_fft(10:floor(end/2))) )
%  imag( (T_z_0_fft(10:floor(end/2)))-(T_z_fft(10:floor(end/2))) )
%
%  sum(sqrt(real( (T_z_0_fft(10:floor(end/2)))-(T_z_fft(10:floor(end/2))) ).^2 + ...
%      imag( (T_z_0_fft(10:floor(end/2)))-(T_z_fft(10:floor(end/2))) ).^2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

err_fft =  sum(abs(log(T_z_0_fft(1:floor(end/2)))- abs(log(T_z_fft(1:floor(end/2)))) ))*1e-2; %+ ...
%  sum(abs( angle(T_z_0_fft(10:floor(end/2))) -
%  angle(T_z_fft(10:floor(end/2))) ));8
%err_E = sum((abs(Q_I)).^2)*1e-4;
err_T_z = sum((T_z(k_min:k)-T_z_0(k_min:k)).^2);
%err_Q = abs(var(Q_I(k_min:k)) - (var(Q_D(k_min:k))));
err_u =  abs(var(u)-1e0)*1e6; %sum((mean(u)).^2)*1e3 +
%err_tau = (sqrt(1/tau_I))*1e2;
err =  err_T_z + err_u + err_fft; % + err_E;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure(1);
% subplot(4,1,1:2); hold off;
% plot(T_z_0(5000:5500),'r'); axis manual; hold on;
% plot(T_z_0(5000:5500) + 0.05*randn(501, 1),'b');

%Constraits
%  cons = ((K_D < K_D_min)  || (K_D > K_D_max)) ||...
%      ((tau_D < tau_D_min)  || (tau_D > tau_D_max))||...
%      ((K_I < K_I_min)  || (K_I > K_I_max))||...
%      ((tau_I < tau_I_min)  || (tau_I > tau_I_max))||...
%      ((a_I < a_I_min)  || (a_I > a_I_max))||...
%      ((b_I < b_I_min)  || (b_I > b_I_max));

cons = (d_D>=d_I);% || (tau_D>tau_I);% || (var(u)>0.5e0))

if cons
    err = 1e1000;
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
        global ts k_z c_z % d_D d_I
        %global K_D_min K_D_max tau_D_min tau_D_max K_I_min ...
        %        K_I_max tau_I_min tau_I_max a_I_min a_I_max b_I_min b_I_max ...
        %        u_1_min u_2_min u_3_min d_D_min d_I_min
        global T_o_data T_z_data u_data Q_P_0 T_0

        K_D = optimValues.bestx(1);
        tau_D = optimValues.bestx(2);
        K_I = optimValues.bestx(3);
        tau_I = optimValues.bestx(4);
        a_I = optimValues.bestx(5);
        b_I = optimValues.bestx(6);
        u_1 = optimValues.bestx(7);
        u_2 = optimValues.bestx(8);
        u_3 = optimValues.bestx(9);
        d_D = round(optimValues.bestx(10));
        d_I = round(optimValues.bestx(11));

        k_0 = max([d_I,d_D])+2;
        u = 0*ones(k_0,1);
        Q_I = 0*ones(k_0,1);
        Q_D = 0*ones(k_0,1);
        K_O = 0*ones(k_0,1)>1;
        T_z = T_z_data(1:k_0);
        k_min = 1000;
        sum_T_z = 0*ones(k_0,1);
        %u = (u_1.*u_data(:,1) + u_2.*u_data(:,2)).*u_data(:,3);
        %u = u_data;
        %         T_0 = 295*ones(1,length(u_data(:,1)));
        %         if 1
        %             T_0(3566:4802) = 294.75;
        %             % %u_0(48023:50330) = 295.5;
        %             T_0(5839:6143) = 296.5;
        %             T_0(6138:15425) = 296;
        %             %u_0(122125:137650) = 296.6;
        %             %T_0(16535:17059) = 295;
        %             %u_0(15425:18125) = 295.5;
        %             T_0(15425:end) = 295;
        %         else
        %             T_0(5805:15431) = 296;
        %         end
        T_z_0 = T_z_data;
        Q_P = u_3*Q_P_0;

        for k = k_0:1:(length(u_data(:,1))-d_I-round(tau_I/ts))

            %sum_T_z(k) = sum_T_z(k-1) + (T_z_data(k-1+d_I)-T_z(k-1));
            %u(k) = a_I*(u_1*(T_z_data(k-1+d_I+round(tau_I/ts))-T_z(k-1)) + u_2*sum_T_z(k-1));
            %max_u = 0.1;
            %u(k-1) = u(k-2) + max(min(u(k-1)-u(k-2),max_u),-max_u);
            %sum_T_z(k) = sum_T_z(k-1) + (u_0(k)-T_z(k-1));
            %u(k-1) = u_3*(u_0(k)-T_z(k-1)) + u_2*sum_T_z(k) + u_1*(T_z_data(k-1+d_I)-T_z(k-1));
            %             if abs(T_0(k)-T_z(k-1))>0.5
            sum_T_z(k) = sum_T_z(k-1) + (T_0(k)-T_z(k-1));
            u(k) = u_1*((T_0(k)-T_z(k-1))+ u_2*sum_T_z(k));
            %u(k-d_I) = u_1*(T_0(k)-T_z(k-1)) + u_2*(T_z_data(k-1)-T_z(k-1));
            %             else
            %                 sum_T_z(k) = 0;
            %                 u(k) = 0;
            %             end
            %K_O(k) = ((T_o_data(k-1) > T_z(k-1)) && (T_0(k-1)>T_z(k-1))) || ((T_o_data(k-1) < T_z(k-1)) && (T_0(k-1)<T_z(k-1)));
            Q_I(k) = Q_I(k-1) + ts*(K_I*atan(u(k-1-d_I))*(b_I + (T_o_data(k-d_I-1)-T_z(k-d_I-1)))-Q_I(k-1))/tau_I; %Energy input
            Q_D(k) = Q_D(k-1) + ts*(K_D*(T_o_data(k-d_D-1)-T_z(k-d_D-1))-Q_D(k-1))/tau_D; %Energy disturbance

            T_z(k) = T_z(k-1) + ts*(Q_I(k) + Q_D(k) + Q_P(k))/c_z;  %Zone temperature
        end

        if 1
            subplot(4,1,1:2); hold on;
            plot(T_z(1:end),'r');
            plot(T_z_data(1:length(T_z)),'b')
            plot(T_0,'k')
            xlim('tight')
            legend('Room temperature model','Room temperature measurements','Location','best')
            subplot(4,1,3); hold off;
            plot(Q_I(1:end)/c_z,'r'); hold on;
            plot(Q_D(1:end)/c_z,'b')
            plot(Q_P(1:end)/c_z,'k')
            xlim('tight')
            legend('Energy controler','Energy outside','Location','best')
            subplot(4,1,4); hold off;
            plot(u(1:end),'b'); hold on;
            legend('Control signal','Location','best')
            ylabel('$u$'); xlabel('Time step')
            xlim('tight')



        end

        segments = [2200,2800; 4600,5200; 7500,8000; 9500,10000];
        len = length(T_z);
        figure(31);
        for i = 1:4
            k_seg = segments(i,1):1:segments(i,2);
            subplot(2, 2, i); hold off;
            plot(k_seg, T_z(k_seg),'r'); hold on;
            plot(k_seg, T_z_data(k_seg),'b');
            xlabel('Time step')
            xlim('tight')
            if i ==2
                legend('Room temperature model','Room temperature measurements','Location','best')
            end
        end


        T = ts;   % Sampling period
        Fs = 1/ts;            % Sampling frequency
        L = (length(u_data(:,1))-d_I-round(tau_I/ts));             % Length of signal
        t = (0:(L-1))*T;        % Time vector
        T_z_0_fft = fft(T_z_0(1:k));
        T_z_fft = fft(T_z(1:k));

        if 1
            figure(42); hold off;
            plot(Fs/L*(0:(L/2-1)),log(abs(T_z_0_fft(1:end/2))),'b'); hold on;
            plot(Fs/L*(0:(L/2-1)), log(abs(T_z_fft(1:end/2))),'r');
            legend('Real','Model')
        end

        k_0 = max([d_I,d_D])+2;
        u = 0*ones(k_0,1);
        Q_I = 0*ones(k_0,1);
        Q_D = 0*ones(k_0,1);
        T_z = T_z_data(1:k_0);
        k_min = 1000;
        sum_T_z = 0*ones(k_0,1);

        %T_0 = 295*ones(1,length(u_data(:,1)));
        %T_0(58205:154311) = 296;

        T_z_0 = T_z_data;
        % figure(11); hold off;
        % plot(T_z_data);  hold on;
        % plot(u_0 + Q_P)

        for k = k_0:1:length(u_data(:,1))-d_I

            u(k) = 0;
            %max_u = 0.1;
            %u(k-1) = u(k-2) + max(min(u(k-1)-u(k-2),max_u),-max_u);
            %sum_T_z(k) = sum_T_z(k-1) + (T_z(k-1)-u_0(k));
            %u(k-1-d_I) = u_1*(T_z(k-1)-u_0(k)) + u_2*sum_T_z(k);

            Q_I(k) = Q_I(k-1) + ts*(K_I*atan(u(k-1-d_I))*(b_I + T_o_data(k-d_I-1)-T_z(k-d_I-1))-Q_I(k-1))/tau_I; %Energy input
            Q_D(k) = Q_D(k-1) + ts*(K_D*(T_o_data(k-d_D-1)-T_z(k-d_D-1))-Q_D(k-1))/tau_D; %Energy disturbance

            T_z(k) = T_z(k-1) + ts*(Q_I(k) + Q_D(k) + Q_P(k))/c_z;  %Zone temperature
        end


        if 1

            figure(41);
            subplot(4,1,1:2); hold off;
            plot(T_o_data); hold on;
            plot(T_z(1:end),'r');
            plot(T_z_data(1:length(T_z)),'b')
            xlim('tight')
            legend('Outside temperature measurements','Room temperature model','Room temperature measurements','Location','best')
            subplot(4,1,3); hold off;
            plot(Q_I(1:end),'r'); hold on;
            plot(Q_D(1:end),'b')
            xlim('tight')
            legend('Energy controler','Energy outside','Location','best')
            subplot(4,1,4); hold off;
            plot(u(1:end),'b'); hold on;
            legend('Control signal','Location','best')
            ylabel('$u$'); xlabel('Time step')
            xlim('tight')
            pause(0.0001)

        end
    case 'done'
        % Finalize your plot here
        hold off;
end

end
