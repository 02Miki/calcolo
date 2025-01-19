%% 1

clc
clear
close all


f = @(x) x^4+4*x;

a = -1;
b = 4;

x0 = 2;
x1 = -1;
k=3;

% devo implementare newton

% 0 = f(xk) + f(xk)'*(xk+1-xk)

% xk+1 = xk - f(xk)/f(xk)'

% approssimo f' con le secanti
% f' = f(xk) - f(xk-1) / xk - xk-1

for n = 1:k
    
    df = (f(x1) - f(x0))/(x1-x0);
    x0 = x1;
    x1 = x1 - f(x1)/df;


end


%% 2

clc
clear
close all


f = @(x) x^3 + sin(x);
df = @(x) 3*x^2+cos(x);

a = -2;
b = 4;
x = 2;

for k = 1:6
    x = x - f(x)/df(x)

end

%% 3
clc
clear
close all


f = @(x) x^4 + 2*x^3 - 4*x-8;

a = 0;
b = 4;

for k = 1:3
    x = (a+b)/2;

    if f(x)*f(a)<0
        b = x;
    else
        a = x;
    end
end



%% 4

clc
clear
close all

f = @(x) exp(-x^2/3);

x = 1;
for k = 1:5
    x = f(x)

end



%% extra gruppo
clc
clear
close all

phi = @(x) -(exp(-x^2)*sin(5*x)+0.3)/10;

x0 = -1;

toll = 10^-3;


for k=1:1000
    
    x = phi(x0);
    
    if abs(x-x0) < toll
        k
        break
    end
    x0=x;

end
















