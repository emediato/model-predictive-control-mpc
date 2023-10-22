function [Sf,bf,Si,bi] = determina_oinf(Af,Gamma,Spsi,bpsi,max_iter,tol)

% [Sf,bf,Si,bi] = determina_oinf(Af,Gamma,Spsi,bpsi,max_iter,tol)
%
% max_iter --> No. m�ximo de itera��es
% tol > 0 --> Toler�ncia no teste de redund�ncia
%            (default = 0)

% Inicializa��o
if nargin < 6, tol = 0; end
r = length(bpsi); % No. restri��es
SpsiGamma = Spsi*Gamma;
S = SpsiGamma; b = bpsi;
flag_redund = 0;
i = 1;

while ( (i <= max_iter) && (flag_redund == 0) )
    disp(['Itera��o ' num2str(i) ' / ' num2str(max_iter)])
    Si{i} = S; bi{i} = b; % Update: 01 Sep 2022
    flag_redund = 1; % Se esse flag n�o for mudado, o algoritmo ser� encerrado 
                     % pois todas as novas restri��es ser�o redundantes
    SGAi = SpsiGamma*(Af^i);
    for j = 1:r % Teste de redund�ncia de cada restri��o
        c = SGAi(j,:)';
        d = bpsi(j);
        % A restri��o � redundante se e somente se t(j) <= tol
        t(j) = teste_redundancia(S,b,c,d);
        if t(j) > tol % Se a restri��o n�o for redundante
            flag_redund = 0; % Pelo menos uma restri��o n�o foi redundante
            % Agrega-se essa restri��o �s anteriores
            S = [S;c']; b = [b;d];
        end
    end
    i = i + 1;
end

% Finaliza��o
if (flag_redund == 0) % O algoritmo foi encerrado sem que houvesse redund�ncia das restri��es
    disp('N�o foi poss�vel caracterizar o MAS')
    Sf = []; bf = [];
else
    Sf = S; bf = b;
end