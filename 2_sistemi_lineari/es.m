%% 1
clc
clear
close all
format long


for n = 4:2:12
    H = hilb(n);
    % somma sulle righe per ottenere il vettore unitario come soluzione
    b = sum(H, 2);
    x = H\b;

    solVera = ones(n, 1);
    errore = abs(x-solVera)/abs(solVera)
end


%% 3

clc
clear
close all

N = 100;

for d = [sqrt(N), N, N*sqrt(N)]
    matrice = spdiags([4*ones(N, 1), -1*ones(N, 1), -1*ones(N, 1), -1*ones(N, 1), -1*ones(N, 1)], [0, 1, -1, d, -d], N, N);
    b = sum(matrice, 2);
    [L, U, P] = lu(matrice);
    % La matrice P (matrice di permutazione) serve solo se c'è pivoting,
    % cambia l'ordine delle righe di A
    y = L\(P*b);
    x = U\y;

end







