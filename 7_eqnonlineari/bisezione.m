function [x, k] = bisezione(a, b, f, tol, kMax)

for k= 1:kMax
    x = (a+b)/2;

    if f(a)*f(x) < 0
        b = x;
    elseif f(x)*f(b) < 0
        a = x;
    else
        disp("intervallo non contiene soluzione oppure intervallo contiene due (o più) soluzioni")
        break
    end
    if abs(f(x)) < tol
        break
    end

end

end