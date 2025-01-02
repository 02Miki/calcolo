%% 1

clc
clear
close all


f = @(x) exp(-x^2/3);

[x, iter, xVec] = puntoFisso(f, 10^-5, 10^-5, 100, 1)
disp("-------")

f = @(x) log(x-1)+3;

[x, iter, xVec] = puntoFisso(f, 10^-5, 10^-5, 100, 1.2)

disp("-------")
f = @(x) sin(x);

[x, iter, xVec] = puntoFisso(f, 10^-5, 10^-5, 100, 0.09)

%% 2

clc
clear
close all


% proviamo newtwon

f = @(x) x^2-2;
df = @(x) 2*x;

[x, k] = newton(5, f, df, 10^-5, 10^-5, 100)


%% 1 aggiuntivi
clc
clear
close all


f = @(x) x^2-2;

[x, k] = bisezione(-1, 15, f, 10^-5, 100)


%% 2 aggiuntivi
% il metodo delle secanti è il metodo di newtwon, ma al posto di usare df,
% si usa 

% df = (f(xk)-(fxk_1))/(xk-(xk_1))




