function [x, k] = newton(x0, f, df, tolf, tolr, kMax)

    % 0 = f(xk) + f'(xk)*(x-xk)
    % 
    % x = xk-f(xk)/f'(xk)
    x = x0;
    for k = 1:kMax
        x_1 = x;
        x = x-f(x)/df(x);

        if abs(x-x_1)/abs(x) <=tolr || abs(f(x)) < tolf
            break
        end
    end

end

