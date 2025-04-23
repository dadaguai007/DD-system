function sig = normalization(data)
s = std(data);
m = mean(data);
sig = (data-m)/s;
end