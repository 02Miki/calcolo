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

%% round 2

%% 1
clc
clear
close all

f = @(x) x.^4+4.*x;

df = @(x) 4.*x.^3 + 4;

x0 = 2;

for k=1:6
    x0 = x0-f(x0)/df(x0)

end

%% 2
clc
clear
close all


phi = @(x) 1/(x.^2 + sin(x) + 1);

x0 = 0;
errori = [];
for k = 1:100

    x = phi(x0);
    errori(k) = abs(x-x0);
    
    if abs(x0-x) < 10^-7
        k
        x0 = x;
        break
    end
    x0 = x
end

errori(2:end)./errori(1:end-1)

% approssimo la derivata di phi con le derivate centrate
h = 0.0001;
abs((phi(x+h) - phi(x-h))/(2*h))

% != 0 -> p = 1


%% 3

clc
clear
close all

f = @(x) 1/3 .*x.^3 + 3.*x + 2.*sin(x) + pi;

a = -1;
b = 0;
for i = 1:12
    x = (a+b)/2
    if (f(a)*f(x)) <0
        b =x;
    else
        a = x;
    end

end


%% 4
clc
clear
close all

f = @(x) x.^4 + 4*x;
% non serve perché metodo secanti
% df = @(x) 4.*x.^3 + 4;

x=[-2 2]';

for k = 1:5
    df = (f(x(k+1)) - f(x(k)))/(x(k+1) - x(k));

    x(end+1) = x(k+1) - f(x(k+1))/df 
end




%% round 3

%% 1

f = @(x) exp(-x.^2./3);


kMax = 5;

x0 = 1;


for k = 1:kMax
    x0 = f(x0)


end


%% 2

clc
clear

f = @(x) x.^4 + 2.*x.^3 - 4.*x - 8;


a = 0;
b = 4;

for k = 1:6
    c = (a+b)/2;
    f(a)*c
    if f(a)*f(c) < 0
        b = c;
    else
        a = c;

    end



end

%% 3
clc
clear

f = @(x) log(x) - 2;

df = @(x) 1/x;

x0 = 3;


for k = 1:5
    x0 = x0-f(x0)/df(x0)

end


%% 4

clc
clear

f = @(x) x.^3 .*sin(x);
x0 = 2;

df = @(x) 3.*x^2*sin(x) + x.^3*cos(x);

for k = 1:4

    x0 = x0-f(x0)/df(x0)
end













