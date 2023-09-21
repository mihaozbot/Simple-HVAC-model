
M.N = 0.05;
M.ts = 60*10;
%M.k_z = 20;
M.c_z = 72e4;
M.act_P = 1;

load('HVAC_par.mat') % Import model parameters: par

M.K_D = par(1);
M.tau_D = par(2);
M.K_I = par(3);
M.tau_I = par(4);
M.a_I = par(5);
M.b_I = par(6);
M.u_1 = par(7);
M.u_2 = par(8);
M.u_3 = par(9);
M.d_D = round(par(10));
M.d_I = round(par(11));

clear par

load('data_HVAC.mat') % Import outside and room temperature data: T_o_data T_z_data u_data

Q_P_0 = (0.5+0.5*sin((90:1:length(u_data(:,1))+90-1)*(2*pi)/144));
Q_P_week = [15:144:(length(Q_P_0)-144);(144+15):144:length(Q_P_0)]';
for iwHVAC = 1:1:length(Q_P_week(:,1))
    if (sum((T_z_data(Q_P_week(iwHVAC,1):Q_P_week(iwHVAC,2)) - min(T_z_data(Q_P_week(iwHVAC,1):Q_P_week(iwHVAC,2))) - ...
            Q_P_0(Q_P_week(iwHVAC,1):Q_P_week(iwHVAC,2))').^2)> ...
            sum((T_z_data(Q_P_week(iwHVAC,1):Q_P_week(iwHVAC,2)) - min(T_z_data(Q_P_week(iwHVAC,1):Q_P_week(iwHVAC,2)))).^2))
        Q_P_0(Q_P_week(iwHVAC,1):Q_P_week(iwHVAC,2)) = 0;
    else
        Q_P_0(Q_P_week(iwHVAC,1):Q_P_week(iwHVAC,2)) = Q_P_0(Q_P_week(iwHVAC,1):Q_P_week(iwHVAC,2))*...
            min(3,max(T_z_data(Q_P_week(iwHVAC,1):Q_P_week(iwHVAC,2)) - ...
            max(T_z_data(Q_P_week(iwHVAC,1)),T_z_data(Q_P_week(iwHVAC,2)))));
    end
end

Q_P_data = M.u_3*Q_P_0;

clear iwHVAC Q_P_week Q_P_0 

