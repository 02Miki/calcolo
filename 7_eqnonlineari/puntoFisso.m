function [x, k, xvec] = puntoFisso(phi, tolf, tolr, kMax, x0)
    
    x = x0;
    xvec(1) = x; 
    for k = 1:kMax
        x = phi(x);
        xvec(k+1) = x;



        if abs((x - xvec(k)))/abs(x) <= tolr || abs(phi(x)-x) <=tolf
            break
        end

    end

end

