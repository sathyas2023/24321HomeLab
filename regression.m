clc;
clear all
close all;

D_l = 8*0.0254; % m
D_s = 7*0.0254; % m
H_l = 6.5*0.0254; % m
H_s = 6*0.0254; % m

A_wat_top_l = pi/4*D_l^2; % m^2
A_wat_top_s = pi/4*D_s^2; % m^2
A_pot_wall_l = pi*D_l*H_l; % m^2
A_pot_wall_s = pi*D_s*H_s; % m^2

Tinf = 21; % C

m_water = 0.71; % kg
c = 4184; % J/(kg-K)

M = table2array(readtable('data.xlsx'));

% Large Pot
t_l = M(:,1); % s
T_l = M(:,2); % C
hL_l = M(:,3); % W/(m^2-K)
hD_l = M(:,4); % W/(m^2-K)
n_l = length(t_l);

% Small Pot
t_s = M(:,5); % s
T_s = M(:,6); % C
hL_s = M(:,7); % W/(m^2-K)
hD_s = M(:,8); % W/(m^2-K)
n_s = length(t_s);

T_water = @(a,b,t,Ti) Tinf+(Ti-Tinf)*(exp(-a*t)+(b/a)/(Ti-Tinf)*(1-exp(-a*t)));

%% Large Pot
e_min = inf;
q_in = -1;
T_w_correct = [];

for q_i = 1:0.01:10000 % q_i is q_in being tested
    e = 0;
    T_w_calc_i = zeros(1,n_l);
    for j = 1:1:n_l
        a = 1/(m_water*c)*(hL_l(j)*A_pot_wall_l+hD_l(j)*A_wat_top_l);
        b = q_i/(m_water*c);
        T_w_calc = T_water(a,b,t_l(j),T_l(j));
        T_w_calc_i(j) = T_w_calc;
        e = e + (T_w_calc - T_l(j))^2;
    end
    if e < e_min
        e_min = e;
        q_in = q_i;
        T_w_correct = T_w_calc_i;
    end
end

fprintf('\nLarge Pot: q_in = %.2f W',q_in);
T_w_correct

%% Small Pot
e_min = inf;
q_in = -1;
T_w_correct = [];

for q_i = 1:0.01:10000 % q_i is q_in being tested
    e = 0;
    T_w_calc_i = zeros(1,n_s);
    for j = 1:1:n_l
        a = 1/(m_water*c)*(hL_s(j)*A_pot_wall_s+hD_s(j)*A_wat_top_s);
        b = q_i/(m_water*c);
        T_w_calc = T_water(a,b,t_s(j),T_s(j));
        T_w_calc_i(j) = T_w_calc;
        e = e + (T_w_calc - T_s(j))^2;
    end
    if e < e_min
        e_min = e;
        q_in = q_i;
        T_w_correct = T_w_calc_i;
    end
end

fprintf('\nSmall Pot: q_in = %.2f W\n',q_in);
T_w_correct   
