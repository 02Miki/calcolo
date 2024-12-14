function r = rho(A, B)
    r = max(abs(eig(-inv(B)*A)))

end




