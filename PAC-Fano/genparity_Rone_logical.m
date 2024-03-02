function u = genparity_Rone_logical(inp,state,polyL)
% polyL is PAC conv of 133 --> 1 011 011
% inp为输入 u为输出
new_state = [inp state];
u = mod(sum(new_state(polyL)),2);%卷积计算

