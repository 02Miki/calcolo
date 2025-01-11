%% 2

clc
clear
close all

f = @(t, y) y+4*t-3;
u0 = 0;

u = u0;
T = 7;
h = 0.01;


for t = 0+h:h:T
    u(end+1) = 1/(1-h) * (u(end) + h*(4*t-3))

end

%% 3

clc
clear
close all


f = @(t,y) [y(2); -y(1)];

u0 = [0; 2^-1];

t0 = 0;
T = 3;

[t, u] = ode45(f, [t0, T], u0);

u(end, 1)

%% 4
clc
clear
close all

f = @(t,u) [y(2), y(2) - 2*y(1) - t^3];

h = 0.01;
T = 3;
t0 = 0;

u0 = [1 5];
u = u0;

for t = t0+h:h:T
    % u1 = u(end, 2)/(1-h);
    u2 = (u(end, 2) + h*(-2*u(end, 1)-t^3))/(1-h+2*h^2);
    u1 = u(end, 1) + h*u2;
    u(end+1, :) = [u1, u2]
    
end



