function Lambda = A2_Q2_getLambda(QP)
    if QP == 1
        Lambda = 0.0431;
    elseif QP == 4
        Lambda = 0.2823;
    elseif QP == 7
        Lambda = 3.8730;
    elseif QP == 10
        Lambda = 26.7773;
    else
        Lambda = 0.026 * (2 ^ QP) + 0.113826;
    end
end