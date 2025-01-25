function [U] = laxFrie(U, lambda, a)
% gestisce un passo
% U = medie di cella
% lambda = deltaT/deltaX
%
% l'ultima cella ha come successiva la prima (si chiama condizione periodica)

N = length(U);

% soluzione di tutte le celle a destra
U_dx = [U(2:end); U(1)];

% soluzione celle a sinistra
U_sx = [U(end); U(1:end-1)];


U = 1/2 * (U_dx+U_sx) - 1/2 * lambda * a * (U_dx-U_sx);

end

