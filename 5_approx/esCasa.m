%% 2

clc
clear
close all

f = @(x) 1./(1+x.^2);
a = -5;
b = 5;

n = 13;

spazio = linspace(a, b, n+1)';
p = polyfit(spazio, f(spazio), n);

spazioBello = linspace(a, b);
val = polyval(p,spazioBello);

figure
hold on
plot(spazioBello, val)
plot(spazioBello, f(spazioBello))

legend("approssimazione", "vera")

%% scarabocchi
clc
clear
close all

f = @(x) x.^3 - 2*x^2-1;

f(2)

fzero(f, [2 3])

clc
clear
close all

f = @(x) x.^3 + 5*x+3-sin(pi*x);

f(2)

fzero(f, [0 -1])


clc
clear
close all

f = @(x) x.^3 - x-5;

f(1)

% fzero(f, [0 -1])



