function u = genparity_Rone_logical(inp,state,polyL)
% polyL is PAC conv of 133 --> 1 011 011
% for a given inout inp, output is u 输入位inp 输出为u
new_state = [inp state];
u = mod(sum(new_state(polyL)),2);%卷积计算

