%% 1
clc
clear
close all

f = @(x) exp(-x+2);

[g, f] = quad(f, 3, 8)

% si ok ma hai avuto culo, guardare come implementare simpson con N intervalli

%% 2
clc
clear
close all

a = -3;
b = 5;
f = @(x) log(x.^4 + x.^2 + 1);

q = [];
for N = [6,12]
    h = (b-a)/N
    spazio = linspace(a, b, N+1)
    q(end+1) = h/6*(f(a) + 2*sum(f(spazio(2:end-1))) + 4*sum(f(spazio(1:end-1)+h/2)) + f(b))
end


q(1)-q(2)

% [g, f] = quad(f, -3, 5, 300000000)



%% 3
clc
clear
close all

f = @(x) 3*exp(-x);

spazio = linspace(2,4,51);
trapz(spazio, f(spazio))



%% extra

f = @(x) sin(3.*x)-x.^2
a = 5;
b=9;

spazio = linspace(a,b, 29+1)
h = (b-a)/29;

h/6*(f(a) + 2*sum(f(spazio(2:end-1))) + 4*sum(f(spazio(1:end-1)+h/2)) + f(b))


%% round 2

%% 1
clc
clear
close all

N = 20;

a = 0;
b = 1;
h = (b-a)/N;
f = @(x) exp(-x.^2);

x = linspace(a,b,N+1);

h/6 * (f(x(1)) + 2 * sum(f(x(2:end-1))) + 4 * sum(f(x(1:end-1) + h/2)) + f(x(end)))


%% 3
clc
clear
close all

f = @(x) exp(-x.^2);

a = -1;
b = 1;
intervalli = 12;
x = linspace(a,b,intervalli+1);

trapz(x, f(x))

%% 4
clc
clear
close all

intervalli = 51;
a = 2;
b = 9;
f = @(x) -3.*sin(x) + cos(x);
x = linspace(a, b, intervalli+1);

trapz(x, f(x))

%% round 3


%% 1

f = @(x) exp(-x.^2);

x = linspace(-1, 1, 13);
trapz(x, f(x))



%% 3

intervalli = 90;

x = linspace(4, 8, intervalli+1);

f = @(x) log(x) + x.*(3.*x);

trapz(x, f(x))


%% 4
clear
intervalli = 20;

a = 0;
b = 1;

f = @(x) exp(-x.^2);

h = (b-a)/intervalli;
x = a:h:b


h/6 * (f(a) + 2 * sum(f(x(2:end-1))) + 4 * sum(f(x(1:end-1) + h/2)) +f(end)  )





























