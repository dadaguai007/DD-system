function wout = my_pec(win, alpha, constellation)
wout = win;

Q = decision(wout, constellation);
error = wout - Q;
error_judge = sign(error);
for i = 2:length(wout)-1
    judge = sign(error(i-1)+error(i+1));
    if error_judge(i) ~= judge
        temp = 0;
    else
        temp = error(i-1)+error(i+1);
    end
    wout(i) = wout(i) + alpha*temp/4;
end

Q = decision(wout, constellation);
error = wout - Q;
error_judge = sign(error);
for i = 2:length(wout)-1
    judge = sign(error(i-1)+error(i+1));
    if error_judge(i) ~= judge
        temp = 0;
    else
        temp = error(i-1)+error(i+1);
    end
    wout(i) = wout(i) + alpha*temp/8;
end

Q = decision(wout, constellation);
error = wout - Q;
error_judge = sign(error);
for i = 2:length(wout)-1
    judge = sign(error(i-1)+error(i+1));
    if error_judge(i) ~= judge
        temp = 0;
    else
        temp = error(i-1)+error(i+1);
    end
    wout(i) = wout(i) + alpha*temp/8;
end

end