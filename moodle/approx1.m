%% 1
clc
clear
close all

x = [0 10 19];
y = [3 7 8];

polyval(polyfit(x, y, 2), exp(0.6))

%% 2
clc
clear
close all

f = @(x) exp(-1./x);
a = 1;
b = 2;

n = 6;
k = 0:n;

% valore = (pi*(2.*k+1))./2*(n+1);
% x = -cos(valore)
% numeratore = pi*(2*k+1);
% denominatore = 2*(n+1);
% 
% numeratore/denominatore
% x = -cos(numeratore/denominatore);

x = -cos((pi*(2*k+1))/(2*(n+1)));

z = (a-b)/2.*x+(a+b)/2;

polyfit(z, f(z), 5)

% il primo ha il grado più alto, -> - 0.5247

%% 3
clc
clear
close all

a = 0;
b = 1;

f = @(x) x.*exp(x);

finto = linspace(a, b, 7);
vero = linspace(a,b);

fit = polyfit(finto, f(finto), 6);
val = polyval(fit, vero)


% plot(vero, f(vero), finto, val)
% legend("vero", "finto")

f(0.3) - polyval(fit, 0.3)
f(0.5) - polyval(fit, 0.5)



%% 4
clc
clear
close all

f = @(x) 1./(x.^4+1);
a = -4;
b = 6;
brutto = linspace(a, b, 7);
bello = linspace(a, b);

p = polyfit(brutto, f(brutto), 6);
val = polyval(p, bello)

plot(bello, f(bello), bello, val)

% bastava la teoria --> 5




%% 5
clc
clear
close all



f = @(x) cos(x);
a = 0;
b = pi;
bello = linspace(a, b);
brutto = linspace(a, b, 5);


polyval(polyfit(brutto, f(brutto), 4), pi/7)

%% round 2

%% 1
clc
clear
close all

f = @(x) x.^2+sin(x.^3);

a = 0;
b = 1;
x = linspace(a,b, 5);

fit = polyfit(x, f(x), 4);

spazio = linspace(0,1);
val = polyval(fit,spazio);


abs(f(0.4)- polyval(fit, 0.4))

abs(f(0.5)- polyval(fit, 0.5))


%% 2

clc
clear
close all

f = @(x) cos(x);

a = 0;
b = 2*pi;

x = linspace(a,b,3);

fit = polyfit(x, f(x), 2);

val = polyval(fit, pi/8)


%% 3

clc
clear
close all

x = [3 7 13];
y = [1 4 10];

fit = polyfit(x, y, 2);
val = polyval(fit, sinh(1.2))


%% 4
clc
clear
close all

t = [0 9 13];
v = [20 31 50];

fit = polyfit(t, v, 2);

val = polyval(fit, 10)

%% 5
clc
clear
close all


f = @(x) sin(x);

a = 0;
b = 2*pi;

x = linspace(a,b,6);

fit = polyfit(x, f(x), 5);
val = polyval(fit, 2*pi/9)


























