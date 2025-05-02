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

%% round 2

%% 1
clc
clear
close all

f = @(t,y) -y(1) + t^2;

[t,y] = ode45(f, [0,7], 0)


%% 2

clc
clear
close all

h = 0.1;

f = @(t,y) -y.^2+t.^2;
tf = 10;

y = pi;
% deve essere -h perché eulero esplicito calcola la roba del passo
% successivo con quelli del passo attuale, quindi usando tf-h come passo
% "attuale", calcolo y di tf
for k = 0:h:tf-h
    
    y = y + h*f(k, y)
end


%% 3

clc
clear
close all

h = 0.3;

f = @(x,y) (sin(x) - 5.*y.^2)/3;

y = 0;
for k = 0.3:h:1.2-h
    k
    y(end+1) = y(end) + h*f(k,y(end))
    
end


%% 4 HEUN ????

%% round 2

%% 3

clc
clear


f = @(t,y) y^2 - t;

y0 = [0];


[t,y] = ode45(f, [0 20], y0)



%% 4
clc
clear

f = @(t,y) 2.*y - 6.*t + 7;

y = 1;

h = 0.1;
t = 0;

for k = 0:h:7*h-h
    y = y + h * f(k, y)
    t = t+1

end
































