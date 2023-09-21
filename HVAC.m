function [T_z_k, Q_I_k, Q_D_k] = HVAC(T_z_k_1, T_z_k_1_d_I, T_z_k_1_d_D, U_k_1_d, T_o_k_1_d_I, T_o_k_1_d_D,  Q_I_k_1, Q_D_k_1, Q_P_k,  M)

    Q_I_k = Q_I_k_1 + M.ts*(M.K_I*atan(U_k_1_d)*(M.b_I + T_o_k_1_d_I-T_z_k_1_d_I)-Q_I_k_1)/M.tau_I; %Energy input
    Q_D_k = Q_D_k_1 + M.ts*(M.K_D*(T_o_k_1_d_D-T_z_k_1_d_D)-Q_D_k_1)/M.tau_D; %Energy disturbance
    T_z_k = T_z_k_1 + M.ts*(Q_I_k + Q_D_k + M.act_P*Q_P_k)/M.c_z ;  %Zone temperature

end
