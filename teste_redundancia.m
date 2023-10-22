function t = teste_redundancia(S,b,c,d)

% A restri��o c'*x <= d ser� redundante se e somente se t <= 0

Sa = [S;c'];
ba = [b;d+1]; % Para limitar superiormente o custo

x = linprog(-c,Sa,ba);
t = (c'*x - d);