function res=BuidSin_CosInput(input_X,a,type)
res = [];
if strcmpi(type,'Sin')
    for i=1:length(input_X)
        res=[res sin(a*input_X(i))];
    end

elseif strcmpi(type,'Cos')
    for i=1:length(input_X)
        res=[res cos(a*input_X(i))];
    end
end
