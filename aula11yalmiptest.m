%{
cd YALMIPfolderShouldbeHere
urlwrite('https://github.com/yalmip/yalmip/archive/master.zip','yalmip.zip');
unzip('yalmip.zip','yalmip')
addpath(genpath([pwd filesep 'yalmip']));
savepath
from: https://yalmip.github.io/tutorial/installation/

%}

% Does YALMIP work at all? If not, we might not even be able to create a variable
x = sdpvar(1)

% Can any solver be called?
optimize(x>= 0, x,sdpsettings('debug',1))

% Problems with a specific solver?
optimize(x>= 0, x,sdpsettings('debug',1,'solver','thissolver'))
