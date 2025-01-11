%% 1
clc
clear
close all


f = @(t, y) -y^2+t;

[t,y] = ode45(f, [0, 2], pi)


%% 3

clc
clear
close all


f = @(t, y) -3/2*y+cos(t)/2;

y = 0;
h = 0.3;

finale = 1.5;

% 1.5-0.3 solo così che l'ultimo valore di y sia y(1.5), perché ci fosse
% anche 1.5 come valore di k, allora l'ultimo valore sarebbe y(1.8)

for k = 0.3:h:(finale-h)
    y(end+1) = y(end) + h*f(k, y(end))

end

%% 4

clc
clear
close all



f = @(t, y) y-8*t;

y = 3;
h = 0.01;
N = 7;

finale = h*N;

for k = 0+h:h:finale
    y(end+1) = 1/(1-h*(-8*k))*y(end)

end







