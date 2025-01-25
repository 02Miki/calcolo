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



