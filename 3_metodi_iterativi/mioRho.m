
function r = mioRho(A, B)

    r = max(abs(eig(A\B)));
    if r < 1
        disp("Converge")
    else
        disp("Non Converge")
    end
end




