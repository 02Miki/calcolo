clc
clear
close all


A =  [1 -1 1; 1 1 1; 1 3 9; 1 2 4];
y = [2;3;-5;7];


ATA = A'*A;
Ay = A'*y;

x = ATA\Ay
sol = [11/2; 9/4; -7/4]

%%


clc
clear
close all


A =  [1 0; 1 1/2; 1 1; 1 2]
y = [0; 3/2; 3/2; 5/2];


ATA = A'*A;
Ay = A'*y;

x = ATA\Ay


%%

clc
clear
close all

x = linspace(-1, 1, 10000);
y = @(x) abs(x);

figure
hold on
plot(x, y(x));

plot(-1/6, (-1/6+1/2)/2*y(-1/6), "or")
plot(-1/2, (-1/6+1/2)/2*y(-1/2), "or")
plot(-1/6*ones(100,1), x)
plot(-1/2*ones(100,1), x)
plot(-(-1/6+1/2)/2, 0, "or")


%%

clc
clear
close all

f = @(x) -3.*cos(x).*(1-3.*sin(x).^2)+2.*(1-2.*sin(x).^2)

x = linspace(0, pi);

plot(x, f(x))

%%
clc
clear
close all

h = pi/2;

y = @(x) cos(x)^3+sin(x)^2

f = @(x) h/6 * (y(x) + 4*y(x+h/2) + y(x+h))

f(0) + f(pi/2)



