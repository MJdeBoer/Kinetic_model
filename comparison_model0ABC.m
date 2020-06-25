% Numerical calculation of the steaddy-state concentratation of states Xi for model 0, A, B and C. 
%
% Marijn de Boer, University of Groningen
% 23/06/2020
%% clear workspace
close all
clear all

%% Calculate transport rates
% Random model parameters
k = 10.^(-3*ones(10,1) + 6.*rand(10,1));     % rate constants k(i) corresponds to k_i
L = 10.^(-2 + 4*rand(1));                    % non-cognate substrate 
l = 10.^(-2 + 4*rand(1));                    % cognate substrate
b = 10.^(-2 + 4*rand(1));                    % SBP
r = 10.^(-2 + 4*rand(1));                    % translocator

% check input
for i = 1:length(k)
    assert(k(i) >= 0, 'k must be positive')
end
assert(L >= 0, 'L must be positive')
assert(l >= 0, 'l must be positive')
assert(b >= 0, 'b must be positive')
assert(r >= 0, 'r must be positive')

% Initialize root finding with fsolve. x0(i) corresponds to the concentration of state Xi
x0 = [b; min([b,l]); min([b,r,l]); min([b,r,l]); r; l; L; min([b,L]); min([b,L,r])]';

% Options for fsolve
options = optimoptions(@fsolve,'TolFun',1e-14,'TolX',1e-14,'MaxFunEvals',4000,'MaxIter',1000);

% Steady-state concentration of the Xi's of model 0:
fun_model0 = @(x) transport_model0(k,x,l,r,b);
x_model0 = fsolve(fun_model0,x0(1:6),options)   % x_model0(i) is X_i^0 as defined in maintext.
check_solution(x_model0,'model 0');

% Steady-state concentration of the Xi's of model A:
fun_modelA = @(x) transport_modelA(k,x,l,r,b,L);
x_modelA = fsolve(fun_modelA,x0(1:8),options)   % x_modelA(i) is X_i^A as defined in maintext.
check_solution(x_modelA,'model A');

% Steady-state concentration of the Xi's of model B:
fun_modelB = @(x) transport_modelB(k,x,l,r,b,L);
x_modelB = fsolve(fun_modelB,x0,options)        % x_modelB(i) is X_i^B as defined in maintext.
check_solution(x_modelB,'model B');

% Steady-state concentration of the Xi's of model C:
fun_modelC = @(x) transport_modelC(k,x,l,r,b,L);
x_modelC = fsolve(fun_modelC,x0,options)        % x_modelC(i) is X_i^C as defined in maintext.
check_solution(x_modelC,'model C');

% Steady state transport rate as defined by Eq. 1 in the maintext
v = (k(6)/r).*[x_model0(4), x_modelA(4), x_modelB(4), x_modelC(4)]

% Relative transport rate as defined by Eq. 2 in the maintext
j = [x_modelA(4), x_modelB(4), x_modelC(4)] ./x_model0(4) 


%% functions
function F = transport_model0(k,x,l,r,b)

F(1) = x(2) + x(3) + x(4) + x(6) - l;
F(2) = x(3) + x(4) + x(5) - r;
F(3) = x(1) + x(2) + x(3) + x(4) - b;
F(4) = k(1)*x(1)*x(6) + k(4)*x(3) - (k(2)+k(3)*x(5))*x(2);
F(5) = k(3)*x(5)*x(2) - (k(4)+k(5))*x(3);
F(6) = k(5)*x(3) - k(6)*x(4);

end

function F = transport_modelA(k,x,l,r,b,L)

F(1) = x(2) + x(3) + x(4) + x(6) - l;
F(2) = x(3) + x(4) + x(5) - r;
F(3) = x(1) + x(2) + x(3) + x(4) + x(8) - b;
F(4) = k(1)*x(1)*x(6) + k(4)*x(3) - (k(2)+k(3)*x(5))*x(2);
F(5) = k(3)*x(5)*x(2) - (k(4)+k(5))*x(3);
F(6) = k(5)*x(3) - k(6)*x(4);
F(7) = k(7)*x(7)*x(1) - k(8)*x(8);
F(8) = x(7) + x(8) - L; 

end

function F = transport_modelB(k,x,l,r,b,L)

F(1) = x(2) + x(3) + x(4) + x(6) - l;
F(2) = x(3) + x(4) + x(5) + x(9) - r;
F(3) = x(1) + x(2) + x(3) + x(4) + x(8) + x(9) - b;
F(4) = k(1)*x(1)*x(6) + k(4)*x(3) - (k(2)+k(3)*x(5))*x(2);
F(5) = k(3)*x(5)*x(2) - (k(4)+k(5))*x(3);
F(6) = k(5)*x(3) - k(6)*x(4);
F(7) = k(7)*x(7)*x(1) + k(10)*x(9) - (k(8)+k(9)*x(5))*x(8);
F(8) = x(7) + x(8) + x(9) - L; 
F(9) = k(9)*x(5)*x(8) - k(10)*x(9);

end

function F = transport_modelC(k,x,l,r,b,L)

F(1) = x(2) + x(3) + x(4) + x(6) - l;
F(2) = x(3) + x(4) + x(5) + x(9) - r;
F(3) = x(1) + x(2) + x(3) + x(4) + x(8) + x(9) - b;
F(4) = k(1)*x(1)*x(6) + k(4)*x(3) - (k(2)+k(3)*x(5))*x(2);
F(5) = k(3)*x(5)*x(2) - (k(4)+k(5))*x(3);
F(6) = k(5)*x(3) - k(6)*x(4);
F(7) = k(7)*x(7)*x(1) + k(10)*x(9) - k(9)*x(5)*x(8);
F(8) = x(7) + x(8) + x(9) - L; 
F(9) = k(9)*x(5)*x(8) - k(10)*x(9);

end

function check_solution(x, y)
    if sum(x > 0) ~= length(x)
        s = ['Attention!, the soltion of '  y  ' contains negative numbers.']; 
        disp(s)
        disp('Try to change the initial conditions x0.')
    end
end
