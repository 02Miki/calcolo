%% 1
clc
clear
close all

h = 0.2;

f = @(x) x.^3 + 1;

der = @(x) (f(x)-f(x-h))/h

der(1.6)

%% 2
clc
clear
close all

h = 0.7;

f = @(x) x^2+3*x-1;

der = @(x) (f(x+h)-f(x))/h

der(1.4)


%% 3

clc
clear
close all

h = 0.6;

f = @(x) x^4+2*x^3-5;

der = @(x) (f(x+h)-f(x-h))/(2*h)

der(1.7)


%% round 2

%% 1
clc
clear
close all

f = @(x) x^3+1;

h = 0.1;

x = 1.3;

(f(x+h)-f(x))/h

%% 2
clc
clear
close all

f = @(x) x^4 + 2*x^3-5;

h = 0.6;
x = 2.5;

(f(x+h) - f(x-h))/(2*h)

%% 3
clc
clear
close all

f = @(x) x^3-x;

h = 0.7;
x = 2.6;

(f(x)-f(x-h))/h













