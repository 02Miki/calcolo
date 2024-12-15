% Jacobi è un metodo iterativo che approssima la soluzione di un 
% sistema lineare Ax=b creando una successione x^(k) che sotto 
% opportune ipotesi tende al vettore soluzone x.
% Gli input per il metodo saranno quindi A,b e la condizione iniziale
% x = x^(0). Come tutti i metodi iterativi si fermarà ad un criterio di 
% stop che coinvolge una tolleranza toll oppure per aver raggiunto 
% un numero massimo di iterazioni itermax

% Il criterio di stop che implementiamo è basato sul residuo relativo: 
% ||Ax^(k)-b||/||b|| < toll  (chiamiamo tale quantità res_rel)

% L'ouput sara x = x^(k), ovvero l'approssimazione che si ha al passo k
% quando il criterio di stop è soddisfatto oppure quando si raggiunge
% itermax. Buona prassi nei metodi iterativi è quella di farsi dare come
% output anche il numero di iterazioni k e il residuo relativo res_rel, 
% ovvero ||Ax^(k)-b||/||b||.

% Riassumendo: OUTPUTS: x, res_rel,k
%              INPUTS: A,b,x,toll, itermax (x all'inizio come input sarà
%              x^(0) che verrà di volta in volta sovrascritta con x =
%              x^(1), poi x^(2) etc...

function [x, res_rel, k] = myJacobi(A,b,x,toll,itermax)

% Estraiamo la dimensione della matrice n
n = size(A,1); % attenzione A deve essere un numero

% Jacobi A = D+E+F ---> posto E+F = N, (NOTA N=A-D) scriviamo il sistema 
% come (D+N)x = b ----> Dx = -N x +b ----> 
% iterativamente Dx^(k+1)=-Nx^(k)+b
% Nel termine -Nx^(k)+b è tutto noto e quindi per trovare x^(k+1) possiamo 
% risolvere un sistema del tipo D x^(k+1) =c, dove c = -Nx^(k)+b

% Estraiamo la matrice diagonale D da A memorizzandola in formato sparso
% con D = spdiags(diag(A),0,n,n); diag(A) estrae la diagonale da A e 
% spdiags memorizza in formato sparso diag(A) nella diagonale principale
% (indice 0). Riscordarsi di passare a spdiags le dimensioni n,n
D = spdiags(diag(A),0,n,n);
% Definiamo la matrice N = A-D 
N = A-D;

% Partiamo con il while: fino a che res_rel > toll && k <= itermax
% Ma PRIMA Inizializzazione: Inizializziamo l'indice di iterazione k ad 1 
% e res_rel in modo che sia "grande"
% mettendo il residuo relativo "grande" siamo sicuri che non 
% si fermi alla prima iterazione
k = 1;
res_rel = 9999;
while res_rel>toll && k<=itermax  %entrambe condizioni soddisfatte (&&)
    % se una delle condizioni non è soddisfatta il ciclo si ferma
      
        % calcolo x ---> risolvo D x^(k+1) =c, dove c = -Nx^(k)+b
        x = D\(-N*x+b);

        % calcolo residuo res_rel: è un residuo relativo
        res_rel = norm(A*x-b)/norm(b); 
        % Attensione la norma di b possiamo 
        % calcolarla una volta sola fuori dal while. 
        % Lo stesso vale per l'inversa di D, vedasi myJacobi_efficient

        % update di k
        k = k +1;
end

end