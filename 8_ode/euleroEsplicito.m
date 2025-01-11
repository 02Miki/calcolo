function [t,u] = euleroEsplicito(f,h, t0, T, u0)

    u = u0;
    
    t = t0:h:T;
    % -1 perché altrimenti avrei x tempi, ma x+1 soluzioni, perché la sol
    % originale la ho già
    for k = t(1:end-1)
        u(end+1) = u(end) + h*f(k, u(end));
    end
end

