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

%%

clc
clear
close all

a = 1/4;
b = 3/4;
dx = 0.5;

f = @(x) x.*(x <=1/2) + 1/2.*(x>1/2);
tutto = integral(f, a,b)/dx

meta = integral(f, 1/4, 1/2)
altrameta = integral(f, 1/2, 3/4)


%% 
clc
clear
close all

t0 = 0;
tf = 1;

x0 = -3/4-2;
xf = 3/4+2;

a = 4;
lambda_a = 3/4;


f = @(x) abs(x).*(abs(x) <=1/2 ) + 1/2.*(abs(x)>1/2);
spazio = linspace(x0,xf);

plot(spazio, f(spazio))
lambda = lambda_a/a;
dx = 0.5;

dt = dx*lambda;

x = x0:dx:xf;

U = [];
for k=1:length(x)-1

    U(end+1,1) = 1/dx * integral(f, x(k), x(k+1))
end

nt = 2;

for t=1:nt
    U_sx = [U(end); U(1:end-1)];
    U_dx = [U(2:end); U(1)];


    U = 1/2 * (U_sx + U_dx) - 1/2 * lambda*a*(U_dx - U_sx)

end


%% 
clc
clear
close all

M = 1000;
f = @(x) cos(cos(x)-sin(x));
x_i = 0;
x_f = pi;
u0 = 1;
nL = -0.5;
h = pi/M;
muf = @(x) exp(sin(x));
mux = h/2:h:pi-h/2;
x = 0:h:pi;

d = [muf(mux(1:end-1))+muf(mux(2:end)) , muf(mux(end))];
cl = [-muf(mux(2:end)), 0];
cu = [0, -muf(mux(2:end))];
A = spdiags([cl' d' cu'],[-1 0 1], M, M)/h^2;
b_f = f(x(2:end))';
b_f(end) = b_f(end)/2;
b_c = zeros(M,1);
b_c(1) = muf(mux(1))./h^2 *u0;
b_c(end) = nL/h;
b = b_f+b_c;

u = A\b;
u =[u0;u];
max(u)
