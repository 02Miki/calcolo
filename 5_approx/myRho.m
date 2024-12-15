% La matrice di iterazione del metodo di Jacobi/Gauss_Siedel è della 
% forma: B = -M^{-1}N
% Quindi input: M,N, output rho

function [rho] = myRho(M,N)

% La matrice di iterazione del metodo di Jacobi/Gauss_Siedel è della 
% forma: B = -M^{-1}N
% Per calcolare B = -M^{-1}N in modo efficiente, anziché calcolare 
% l'inversa della matrice M, risolviamo un sistema del tipo: B = - M \ N;

% B = - inv(M)*N; % inv non è operazione stabile
 B = -M\N;

% che coincide con il risolvere n (n = numero colonne di N) sistemi
% lineari del tipo:
% B(:, i) = - M \ N(:, i) per ogni i =1,...,n


% Il metodo iterativo converge se e solo se il raggio spettrale di 
% B è < 1.

% Calcoliamo il raggio spettrale di B: rho = max(abs(eig(B)));

% rho = max(abs(eig(B))); % eig calcola tutti gli autovalori di 
% B e risulta più costosa di abs(eigs(B, k, 'LM'))
% il quale calcola soltanto i k autovalori
% di B più grandi in modulo (LM = largest in magnitude)

rho = abs(eigs(B,1,'LM'));

end