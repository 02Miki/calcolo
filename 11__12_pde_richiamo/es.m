%% 1 differenze finite
% errore: h dovevi mettere N-1 e non N

clc
clear
close all


mu = 1;

BC = [1 -1];

a = 0;
b = 1;

N =100;

% h deve essere N-1 perché N è il numero di nodi, non il numero di intervalli
h = (b-a)/(N-1)

x = linspace(a,b,N)';

uni = ones(N,1);
d = 2*mu/h^2*uni;
d(1) = d(1)/2;
d(end) = d(end)/2;

% un elemento in più lasciato apposta per spdiags
d1 = mu/h^2 * -uni;

A_mio = spdiags([d, d1, d1], [0 1 -1], N, N);

gamma_f = @(x) (1-2.*x).^2-2; 
A_richiamo = spdiags([gamma_f(x)], 0, N, N);

% dobbiamo dividere per due perché abbiamo diviso per due A per neumann
A_richiamo(1) = A_richiamo(1)/2;
A_richiamo(end) = A_richiamo(end)/2;


b = zeros(N,1);
b(1) = b(1) - BC(1)/h;
b(end) = b(end) + BC(end)/h;

u = (A_mio+A_richiamo)\b;

figure
hold on

plot(x,u, "r")

u_esatta = @(x) exp(x-x.^2);

plot(x, u_esatta(x))


%% elementi finiti
% errori:
% non usato xmedi per dgamma
% 
close all

gamma_f = @(x) (1-2.*x).^2-2; 
mu = 1;

BC = [1 -1];

a = 0;
b = 1;

N =100;


% h deve essere N-1 perché N è il numero di nodi, non il numero di intervalli
h = (b-a)/(N-1);
x = linspace(a,b,N)';
uni = ones(N,1);
d = 2*mu/h*uni;
d(1) = d(1)/2;
d(end) = d(end)/2;

d1 = -mu/h*uni;


A = spdiags([d, d1, d1], [0 1 -1], N, N);

b = zeros(N,1) + 2*h;

b(1) = h/2 - BC(1);
b(end) = h/2 + BC(2);
xMedi = (x(2:end) + x(1:end-1))/2;

dgamma = 1/3 * (gamma_f(xMedi(1:end-1)) + gamma_f(xMedi(2:end)))*h;
dgamma = [gamma_f(xMedi(1))/3*h; dgamma; gamma_f(xMedi(end))*h/3];


dgamma1 = 1/6 * gamma_f(xMedi(1:end)) *h;

A_gamma = spdiags([dgamma, [0;dgamma1], [dgamma1;0]], [0 1 -1], N, N);

b = zeros(N, 1);

b(1) = b(1) * h/2 - BC(1);
b(end) = b(end) * h/2 + BC(2);



u = (A+A_gamma)\b;
figure
hold on

plot(x,u, "r")

u_esatta = @(x) exp(x-x.^2);

plot(x, u_esatta(x))


%% 2
% errore: !! h = b-a, non a-b

close all

a = 0;
b = 1;
% metto 100 nodi
N = 100;
BC = [log(5) -log(5)];
mu = 1;

fgamma = @(x) (1-2.*x).^2.*log(5)^2 - 2.*log(5);

% N-1 perché con 100 nodi ho 99 intervalli
h = (b-a)/(N-1);
x = linspace(a, b, N)';
uni = ones(N,1);

% d = mu/h + mu/h
d = 2*mu/h.*uni;
d(1) = d(1)/2;
d(end) = d(end)/2;

% c'è un elemento in più per spdiags
d1 = -mu/h.*uni;

A = spdiags([d, d1, d1], [0 1 -1], N, N);

xMedi = (x(1:end-1) + x(2:end))/2;

dgamma = 1/3 * (fgamma(xMedi(1:end-1)) + fgamma(xMedi(2:end)))*h;
dgamma = [h/6*fgamma(xMedi(1)); dgamma; h/6*fgamma(xMedi(end))];
dgamma_1 = h/6 * fgamma(xMedi);


A_gamma = spdiags([dgamma [0; dgamma_1] [dgamma_1; 0]], [0 1 -1], N, N);

b = zeros(N,1);

b(1) = b(1)*h/2 - BC(1);
b(end) = b(end)*h/2 + BC(2);



u = (A+A_gamma)\b;

figure
hold on
plot(x,u)

u_vera = @(x) 5.^(x-x.^2);

plot(x, u_vera(x))

legend("calcolata", "vera")

%% 3 -- videlezione 16/12

clc
clear
close all


a = 0;
b = 1;
T = 10;
BC = [0 1];
rho = 1;
f = @(t) 1+t.^2;
mu = 5;

f_u0 = @(x) x.^2;

Tf = 1;
T0 = 0;

% nodi interni
N = 100;
% metto h variabile
x = [0; sort(rand(N,1)); 1];
h = diff(x);

% non ci sono coefficienti variabili né per matrice rigidezza né quella di
% massa

% abbiamo dirichlet agli estremi, quindi non si considerano i nodi esterni
Ndof = N;

d = mu*(1./h(1:end-1) + 1./h(2:end));

d1 = - mu./h(2:end-1); 

A = spdiags([d, [0; d1], [d1; 0]], [0 1 -1], N, N);

% matrice di massa (quella davanti al tempo) -> rho * (h/3 + h/3)
dM = rho * (h(1:end-1)/3 + h(2:end)/3);

% codiag è rho*h/6
d1M =rho * h(2:end-1)/6;

M = spdiags([dM, [0; d1M], [d1M; 0]], [0 1 -1], N, N);

% condizioni iniziali

% valutazione della funzione al tempo t = 0, per tutte le x
% devo togliere il primo e l'ultimo elemento perché sono esterni e li so già
U0 = f_u0(x(2:end-1));

% termine noto (dipenderà ancora dal tempo, che si lascia come variabile)
% la prima parte è formula quadratura dei trapezi solita

% tra quadrate c'è contributo dirichlet (lo zeros dentro le quadrate è per
% non andare a cambiare nulla nei nodi interni, ma ci serve perché le
% parentesi quadrate devono essere un vettore completo
b = @(t) f(t)/2 * (h(1:end-1) + h(2:end)) + [mu*BC(1)/h(1); zeros(N-2, 1); mu*BC(2)/h(end)];


% sono in questa situazione:

% M*u' + A*u = q
% u' = (q-A*u)*M^-1
% u' = M\(q-A*u)

% per EE deltaT < 2 / (max(abs(eig(M\A)))

dt = 2/max(abs(eigs(M\A)));


% Si fa EI:
% u' = M\(q-A*u)
% u+1 = u + u'
% u+1 = u + h*(q-A*u+1)*M^-1
% u+1*(M + A*h) = u*M + hq

% u+1 * matrice = tn

dt = 0.1;
Nstep = (Tf-T0)/dt;
T(1) = T0;
matrice = M+A*dt;
U = U0;
for k = 1:Nstep
    T(end+1) = T(end) + dt;
    tn = M*U(:,end)+dt*b(T(k+1));

    U(:,end+1) = matrice\tn;

end
% devo aggiungere dirichlet alla fine, che deve essere l'ultima x di tutti
% i vari tempi
U = [U; BC(2) * ones(size(U,2))];
% U -> sulle righe varia x, su colonne varia t






