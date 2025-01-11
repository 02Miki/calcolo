%% 2 e 3

clc
clear
close all

t0 = 0;
T = 1;

f = @(t,u) u+sin(t);

u0 = 0;

h = 0.1;
[t01,u01] = euleroEsplicito(f, h, t0, T, u0)

h = 0.01;
[t001,u001] = euleroEsplicito(f, h, t0, T, u0)

fvera = @(t) 1/2*(exp(t)-sin(t)-cos(t));

[tMat, uMat] = ode45(f, [t0,T], u0)

figure 
hold on

plot(t01, u01);
plot(t001, u001);
x = linspace(t0, T)
plot(x, fvera(x));
plot(tMat, uMat);


legend("t01", "t001", "vera", "ode45")


%% 4

clc
clear
close all


f = @(t, u) -10*u;
t0 = 0;
T = 10;
u0 = 1;

h = 0.3
[t03,u03] = euleroEsplicito(f, h, t0, T, u0)

h = 0.1
[t01,u01] = euleroEsplicito(f, h, t0, T, u0)


figure 
hold on

plot(t01, u01);
plot(t03, u03);


legend("t01", "t03")

f = @(t, u) -5*u;

h = 0.3
[t03,u03] = euleroEsplicito(f, h, t0, T, u0)

h = 0.1
[t01,u01] = euleroEsplicito(f, h, t0, T, u0)

figure 
hold on

plot(t01, u01);
plot(t03, u03);
legend("t01", "t03")


%% 5 - Robertson
clc
clear
close all


f = @(t,y) [-0.04*y(1) + 10^4*y(2)*y(3); 0.04*y(1) - 10^4*y(2)*y(3) - 3*10^7*y(2)^2; 3*10^7*y(2)^2]

a0 = 0;
a = 40;

u0 = [1; 0; 0]

[tode45, uode45] = ode45(f, [a0, a], u0);

figure
hold on
set(gca,'YScale','log','XScale','log')
plot(tode45, uode45)

%% 5 Van der pol
clc
clear
close all


mu = 2;
f = @(t,y) [y(2); mu*(1-y(1)^2)*y(2)-y(1)];


u0 = [2; 0];


[t, u] = ode15s(f, [0, 30], u0);
figure
hold on
plot(t, u(:, 1))


plot(t, u(:, 2))
xlabel("Tempo")
ylabel("Spazio/Tempo")
legend("Spazio", "Velocità")

%% es aggiuntivi

clc
clear
close all

u0 = [0; 1];

t0 = 0;
T = 25;

h = 0.01;

f = @(t, y) [y(2); -2*y(1) - 3*y(2)]


[t, u] = euleroEsplicitoSistemi(f, h, t0, T, u0)


figure
hold on
plot(t, u(:, 1),'r', 'linewidth', 2)
plot(t, u(:, 2),'b', 'linewidth', 2)
title('EE')
