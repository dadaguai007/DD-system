function y_d=DFE_decision(yn)


% pam4 的判定
if yn > 2
    y_d = 3;
elseif yn > 0
    y_d = 1;
elseif yn > -2
    y_d = -1;
else
    y_d = -3;
end


% Reliability calculation

if abs(yn)>3
    r=1;
else
    r=1-abs(yn-y_d);
end

% IWFE
y_out=yn+r*(y_d-yn);


% 针对不同情况，需要自己写判定程序
end