function [Sf,bf,Si,bi] = determina_oinf(Af,Gamma,Spsi,bpsi,max_iter,tol)

% [Sf,bf,Si,bi] = determina_oinf(Af,Gamma,Spsi,bpsi,max_iter,tol)
%
% max_iter --> No. máximo de iterações
% tol > 0 --> Tolerância no teste de redundância
%            (default = 0)

% Inicialização
if nargin < 6, tol = 0; end
r = length(bpsi); % No. restrições
SpsiGamma = Spsi*Gamma;
S = SpsiGamma; b = bpsi;
flag_redund = 0;
i = 1;

while ( (i <= max_iter) && (flag_redund == 0) )
    disp(['Iteração ' num2str(i) ' / ' num2str(max_iter)])
    Si{i} = S; bi{i} = b; % Update: 01 Sep 2022
    flag_redund = 1; % Se esse flag não for mudado, o algoritmo será encerrado 
                     % pois todas as novas restrições serão redundantes
    SGAi = SpsiGamma*(Af^i);
    for j = 1:r % Teste de redundância de cada restrição
        c = SGAi(j,:)';
        d = bpsi(j);
        % A restrição é redundante se e somente se t(j) <= tol
        t(j) = teste_redundancia(S,b,c,d);
        if t(j) > tol % Se a restrição não for redundante
            flag_redund = 0; % Pelo menos uma restrição não foi redundante
            % Agrega-se essa restrição às anteriores
            S = [S;c']; b = [b;d];
        end
    end
    i = i + 1;
end

% Finalização
if (flag_redund == 0) % O algoritmo foi encerrado sem que houvesse redundância das restrições
    disp('Não foi possível caracterizar o MAS')
    Sf = []; bf = [];
else
    Sf = S; bf = b;
end