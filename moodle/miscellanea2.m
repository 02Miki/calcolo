%% 1

clc
clear
close all


N = 20;

A = zeros(N,N);
for m = 1:N

    for k = 1:N
        A(m,k) = 1 +exp(-abs(m-k));

    end
end



b = zeros(N, 1);
for k=1:2:20
    b(k,1) = 1;
end

x = A\b;
norm(x, inf)

%% 2
clc
clear
close all


f = @(x) log(x.^4 + x.^2 + 1);

N = 6;
N2 = 12;

a = -3;
b = 5;

simp = @(h,x) h/6*(f(x(1)) + 2 * sum(f(x(2:end-1))) + 4*sum(f(x(1:end-1) + h/2)) + f(x(end)));


h = (b-a)/N;
h2 = (b-a)/N2;


x = a:h:b;
x2 = a:h2:b;

abs(simp(h,x) - simp(h2,x2))

%% 3
clc
clear
close all



f = @(x) exp(-1./(x+1));

a = 0;
b = 1;
N = 6;
k = 1:N;
x = -cos((pi*(2*k+1))/(2*(N+1)));

z = (a+b)/2 + (b-a)/2*x;

fit = polyfit(z, f(z), 5)


%% 3

clc
clear
close all


f = @(x) x*log(x^2+3) + sqrt(x^2+4) - 1;
df = @(x) log(x^2+3) + (2*x^2)/(x^2+3)+x/(sqrt(x^2+4));

x = 10;
for k = 1:3
    x = x - f(x)/df(x)

end


%% 4

clc
clear
close all


f = @(t,y) [y(2); y(3); -5*y(3)-8*y(2) - 4*y(1) + t^2+t+1];

[t,y] = ode45(f, [1,6], [3; 4; 5]);

y(end, 2)



%% 5 

clc
clear
close all


a = 0.5;
lambda = 1;
dx = 0.1;
dt = lambda*dx;

t0 = 0;
tf = 1;

x0 = 0;
xf = 1;

f = @(x) sin(2*pi.*x);

x = x0:dx:xf;

for k = 1:length(x)-1
    
    U(k,1) = 1/dx * integral(f, x(k), x(k+1));

end

nt = (tf-t0)/dt;


for k = 1:nt
    U_sx = [U(end); U(1:end-1)];

    U = U - lambda*a * (U - U_sx);

end

max(U)




