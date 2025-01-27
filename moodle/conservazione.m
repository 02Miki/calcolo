%% 1

clc
clear
close all

f = @(x) cos(2*pi.*x);

Tf = 1;
T0 = 0;
lambda = 1;
a = 0.5;
dx = 0.01;

dt = dx;

passiT = (Tf-T0)/dt;

% mi serve un numero di passi intero (in questo caso non serve)

passiT = ceil(passiT);

% ricalcolo dt (dt nuovo <= dt vecchio)

dt = (Tf-T0)/passiT;

% calcolo le uj di partenza

x0 = 0;
xF = 1;
x = x0:dx:xF;
U = [];
for k=1:length(x)-1
    U(k,1) = 1/dx * integral(f, x(k), x(k+1));
end

xMedi = (x(2:end) + x(1:end-1))/2;
figure
for k=1:passiT
    U_dx = [U(2:end); U(1)];
    U_sx = [U(end); U(1:end-1)];

    U = U - 1/2*lambda*a*(U_dx - U_sx) + 1/2 * (lambda*a)^2 * (U_sx + U_dx - 2*U);
    plot(xMedi, U, "bo-")
    pause(0.01)
end

max(U)


%% 2

clc
clear
close all

a = 0.5;

lambda = 1;

dx = 0.1;
dt = lambda*dx;

tf = 1;
f0 = 0;

nt = tf/dt;

x0 = 0;
xF = 1;
x = x0:dx:xF;
f = @(x) sin(2*pi.*x);

for k = 1:length(x)-1
    U(k,1) = 1/dx * integral(f, x(k), x(k+1));
end

xMedi = (x(2:end) + x(1:end-1))/2;
figure
for k = 1:nt
    U_dx = [U(2:end); U(1)];
    U_sx = [U(end);U(1:end-1)];

    U = 1/2 * (U_sx + U_dx) - 1/2 * lambda * a * (U_sx - U_dx);
plot(xMedi, U, "bo-")
    pause(0.01)
end


max(U)


%% 3
clc
clear
close all


a = 0.5;

t0 = 0;
tf = 10;

x0 = 0;
xf = 1;

lambda_a = 0.5;

lambda = lambda_a/a;
 
figure
hold on
for dx = [0.1 0.05 0.02 0.01]
    x = x0:dx:xf;
    dt = lambda*dx;
    
    f = @(x) cos(2*pi.*x);
    nt = (tf-t0)/dt;
    
    for k = 1:length(x)-1
        U(k,1) = 1/dx * integral(f, x(k), x(k+1));
    end
    
    xMedi = (x(2:end) + x(1:end-1))/2;
    
    for k = 1:nt
        U_sx = [U(end);U(1:end-1)];
        U_dx = [U(2:end); U(1)];
    
        U = U - 1/2 * a *lambda * (U_sx-U_dx) + 1/2 * (lambda * a)^2 * (U_dx + U_sx - 2*U);
        
    end
    size(xMedi)
    size(U)
    plot(xMedi, U, "o-")

    clear U

end

x = linspace(0,1);
u_vera = f(x-a*10);

plot(x, u_vera, LineWidth=3)

legend("0.1", "0.05", "0.02", "0.01", "vera")



