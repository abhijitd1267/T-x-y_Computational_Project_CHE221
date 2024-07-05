clc;
clear all;

%%
water_data = readmatrix('water_data.txt');
T = water_data(:,1); % in K
x_j = water_data(:,2); % mole fraction of water
x_i = 1 - x_j; % mole fraction of THF

lambda_ij_ii = 1475.26*4.2; % J/mol
lambda_ji_jj = 1844.79*4.2; % J/mol
V_i = 81.55*10^(-6); % m^3/mol
V_j = 18.07*10^(-6); % m^3/mol

R = 8.314; % J/mol.K
x_set = linspace(0,1,100);

for i = 1:size(T,1)
    WilsonParameter_ij(i,1) = (V_j/V_i)*exp(-lambda_ij_ii/(R*T(i)));
    WilsonParameter_ji(i,1) = (V_i/V_j)*exp(-lambda_ji_jj/(R*T(i)));
end

for i = 1:size(T,1)
    gamma_i(i,1) = exp(-log(x_i(i) + WilsonParameter_ij(i)*x_j(i)) + x_j(i)*(WilsonParameter_ij(i)/(x_i(i)+WilsonParameter_ij(i)*x_j(i)) - WilsonParameter_ji(i)/(x_j(i)+WilsonParameter_ji(i)*x_i(i))));
    gamma_j(i,1) = exp(-log(x_j(i) + WilsonParameter_ji(i)*x_i(i)) - x_i(i)*(WilsonParameter_ij(i)/(x_i(i)+WilsonParameter_ij(i)*x_j(i)) - WilsonParameter_ji(i)/(x_j(i)+WilsonParameter_ji(i)*x_i(i))));
end

% Antoine Equation
A_i = 7.1057;
B_i = 1256.68;
C_i = 232.621;

A_j = 8.07131;
B_j = 1730.63;
C_j = 233.426;

T_new = T - 273; % in deg C

for i = 1:size(T,1)
    Pvp_i(i,1) = (10^(A_i - B_i/(C_i + T_new(i))))*10^5/760; % in Pa
    Pvp_j(i,1) = (10^(A_j - B_j/(C_j + T_new(i))))*10^5/760; % in Pa
end

for i = 1:size(T,1)
    P(i,1) = x_i(i)*gamma_i(i)*Pvp_i(i) + x_j(i)*gamma_j(i)*Pvp_j(i);
end

for i = 1:size(T,1)
    y_i(i,1) = (gamma_i(i)*x_i(i)*Pvp_i(i))/P(i);
    y_j(i,1) = (gamma_j(i)*x_j(i)*Pvp_j(i))/P(i);
end

set_x = [0,0];
set_y = [1,1];

figure(1);
plot(x_i, T);
hold on;
plot(y_i, T);
xlabel('x,y for THF');
ylabel('Temperature (T)');
title('T-x-y curve for THF');
legend('x', 'y');
hold off;

figure(2);
plot(x_j, T);
hold on;
plot(y_j, T);
xlabel('x,y for water');
ylabel('Temperature (T)');
title('T-x-y curve for water');
legend('x', 'y');
hold off;

%% Part 2.
THF_literature_data = readmatrix('THF_literature.txt');
x_THF_lit = THF_literature_data(:,1);
y_THF_lit = THF_literature_data(:,2);

figure(3);
plot(x_i, y_i, '-');
hold on;
plot(x_THF_lit, y_THF_lit, 'og');
plot(x_set, x_set, '--');
xlabel('Mole fraction of THF in liquid phase (x)');
ylabel('Mole fraction of THF in vapor phase (y)');
title('x - y curve for THF');
legend('Simulated value', 'Literature value', 'x = y line');
hold off;

%% Part 3.
T_3 = 50; % in deg C

Pvp_i_50 = (10^(A_i - B_i/(C_i + T_3)))*10^5/760; % in Pa
PVp_j_50 = (10^(A_j - B_j/(C_j + T_3)))*10^5/760; % in Pa

x_i_3 = linspace(0,1,101)';
x_j_3 = 1 - x_i_3;

WilsonParameter_ij_3 = (V_j/V_i)*exp(-lambda_ij_ii/(R*(T_3+273)));
WilsonParameter_ji_3 = (V_i/V_j)*exp(-lambda_ji_jj/(R*(T_3+273)));


for i = 1:size(x_i_3,1)
    gamma_i_3(i,1) = exp(-log(x_i_3(i) + WilsonParameter_ij_3*x_j_3(i)) + x_j_3(i)*(WilsonParameter_ij_3/(x_i_3(i)+WilsonParameter_ij_3*x_j_3(i)) - WilsonParameter_ji_3/(x_j_3(i)+WilsonParameter_ji_3*x_i_3(i))));
    gamma_j_3(i,1) = exp(-log(x_j_3(i) + WilsonParameter_ji_3*x_i_3(i)) - x_i_3(i)*(WilsonParameter_ij_3/(x_i_3(i)+WilsonParameter_ij_3*x_j_3(i)) - WilsonParameter_ji_3/(x_j_3(i)+WilsonParameter_ji_3*x_i_3(i))));
end

for i = 1:size(gamma_i,1)
    P_3(i,1) = x_i_3(i)*gamma_i_3(i)*Pvp_i_50 + x_j_3(i)*gamma_j_3(i)*PVp_j_50;
end

for i = 1:size(x_i_3,1)
    y_i_3(i,1) = gamma_i_3(i)*x_i_3(i)*Pvp_i_50/P_3(i);
    y_j_3(i,1) = 1 - y_i_3(i,1);
end

figure(4);
plot(x_i_3, P_3);
hold on;
plot(y_i_3, P_3);
xlabel('x,y for THF');
ylabel('Pressure (P)');
title('P-x-y curve for THF at T = 50 deg C');
legend('x', 'y');
hold off;

figure(5);
plot(x_i_3, y_i_3);
hold on;
plot(x_set, x_set, '--');
xlabel('Mole fraction of THF in liquid phase (x_1)');
ylabel('Mole fraction of THF in vapor phase (y_1)');
title('x-y diagram for THF at T = 50 deg C');
legend('x-y curve', 'y = x line');
hold off;


