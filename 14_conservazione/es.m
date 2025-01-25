clc
clear
close all

L = 1;
a = 0.5;

% numero intervalli
N = 200;
dx = L/N;

% cour = abs(a) * lambda
% lambda = dT/dX
% lo scelgo in modo tale da garantire che cour < 1 --> dT <= dX/a
% 0.9 numero a caso scelto per provare
dt = dx/abs(a) * 0.9;

Tf = 5;
T0 = 0;

% Molto probabilmente questo valore non è intero
nT = (Tf - T0)/dt;

% lo trasformo in valore intero per eccesso, ma devo verificare che dt <=dx/a
% arrotondando per eccesso, sono sicuro che dt nuovo sia < di dt vecchio
nT = ceil(nT);

% e ricalcolo dt con nT nuovo così che sia tutto coerente
dt = (Tf - T0)/nT;

x = linspace(0,L, N+1)

u0 = @(x) 1.*(0.4 < x & x < 0.6);

U = zeros(N,1);

for k=1:N
    U(k) = 1/dx * integral(u0, x(k), x(k+1));
end

lambda = dt/dx;
% solo per fare il plot
xMedi = (x(2:end) + x(1:end-1))/2;
figure
for k=1:nT
    U = laxFrie(U, lambda, a);
    plot(xMedi, U, "bo-")
    pause(0.001)
end




% LF diffonde (o dissipa), la soluzione infatti tende ad appiattirsi
% LF disperde, la soluzione anticipa o è in ritardo rispetto alla soluzione esatta 

% nel caso di LF, LF tende ad anticipare, cioé è più veloce della soluzione esatta

% Per non avere dissipazione o diffusione, serve cour = 1. In questo caso,
% i tre metodi sono equivalenti


% L'unica differenza degli altri rispetto a LF è che laxwendroff genera
% oscillazioni nei bordi della soluzione (upwind e lf non lo fanno)

% LW è un metodo di ordine superiore











