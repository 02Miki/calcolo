%% 1

clc
clear
close all

x = [0 1 2 3]';
y = [1 2 4 8]';


% a0*1+a1*e^x
A = [ones(length(x),1) exp(x)]

x = A\y

residuo = norm(A*x-y)

% !!!!!! Per far venire il risultato giusto, servirebbe un ^2, ma ho
% chiesto alla prof ed ha detto che è l'esercizio sbagliato, quindi è
% giusto così !!!!!

%% 3
clc
clear
close all

x = [-4 6 16 23];
y = [5 3 4 9];

spline(x, [5 y 6], log(1.3))

%% round 2

%% 1
clc
clear
close all

x = [0 4 12];
y = [6 7 4];

spline(x, [5 y 6], log(0.8))



%% 2

clc
clear
close all

x = [0 1 2 3 4 5]';
y = [1 2 4 8 16 32]';

% a0 + a1 * e^x

A = [ones(length(x),1), exp(x)];

sol = A\y

residuo = norm(A*sol-y)

%% round 3

%% 1


x = [4 6 8 9];
y = [8 6 4 5 9 5];

spline(x,y,exp(0.4))


%% 2

% a0*1 + a1*exp(x) + a2 * exp(2*x)

x = [0 1 2 3]';
y = [1 3 9 27]';

A = [ones(length(x),1), exp(x), exp(2.*x)];

b = A\y

r = norm(A*b-y,2)










