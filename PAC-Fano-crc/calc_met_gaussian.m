function metric = calc_met_gaussian(LLR,u,b)%计算路径度量
Lj = exp(LLR);
if u == 0
    metric = 1 -  log2(1+1/Lj) - b;%0
else
    metric = 1 - log2(1+Lj) - b;%1
end

