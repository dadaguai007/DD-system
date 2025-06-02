function wout = my_pec_advanced(win, alpha, constellation)
ML = length(alpha);

Q = decision(win, constellation);
error = win - Q;
error_judge = sign(error);
wout = win;
for i = 1+ML:length(wout)-ML
    temp = error(i-1:-1:i-ML) + error(i+1:i+ML);
    temp = alpha * temp;
    judge = sign(temp);
    if error_judge(i) == judge
        wout(i) = wout(i) + temp/4;
    end
end

Q = decision(wout, constellation);
error = wout - Q;
error_judge = sign(error);
for i = 1+ML:length(wout)-ML
    temp = error(i-1:-1:i-ML) + error(i+1:i+ML);
    temp = alpha * temp;
    judge = sign(temp);
    if error_judge(i) == judge
        wout(i) = wout(i) + temp/8;
    end
end

Q = decision(wout, constellation);
error = wout - Q;
error_judge = sign(error);
for i = 1+ML:length(wout)-ML
    temp = error(i-1:-1:i-ML) + error(i+1:i+ML);
    temp = alpha * temp;
    judge = sign(temp);
    if error_judge(i) == judge
        wout(i) = wout(i) + temp/8;
    end
end

end