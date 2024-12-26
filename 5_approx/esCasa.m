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

%% 4
clc
clear
close all












