%% Uppgift 1
a = 0.5;
y0 = [1-a;0;0;sqrt((1+a)/(1-a))];
%% d)
N = 10000;
t_end = 100;
h=t_end/N;

%Framåt Euler
y_feuler = fram_euler(y0, N, t_end);
energi_feuler = energi(y_feuler);
size(y_feuler)
q1f = y_feuler(:,1);
q2f = y_feuler(:,2);

figure
plot(0:h:t_end, q1f, 0:h:t_end, q2f, 0:h:t_end, energi_feuler)

%Bakåt Euler
y_beuler = bak_euler(y0, N, t_end);
size(y_beuler)
q1b = y_beuler(:,1);
q2b = y_beuler(:,2);
figure
plot(0:h:t_end, q1b, 0:h:t_end, q2b)

%Symplektisk Euler 
y_seuler = symp_euler(y0, N, t_end);
energi_seuler = energi(y_seuler);
q1s = y_seuler(:,1);
q2s = y_seuler(:,2);
figure
plot(0:h:t_end, q1s, 0:h:t_end, q2s, 0:h:t_end, energi_seuler)
%% Funktioner 
function ydot = kepler(y)
%y = (q1, q2, p1, p2)
q1dot = y(3);
q2dot = y(4);
p1dot = -1*y(1)/(y(1)^2+y(2)^2)^(3/2);
p2dot = -1*y(2)/(y(1)^2+y(2)^2)^(3/2);
ydot = [q1dot;q2dot;p1dot;p2dot];
end

function y_values = fram_euler(y0, N, t_end)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
h = t_end/N;
y = y0;
y_values = [zeros(N+1, 1),zeros(N+1, 1),zeros(N+1, 1),zeros(N+1, 1)];
for i  = 1:N+1
    for k = 1:4
        y_values(i,k) = y(k);
    end
    y = y+h*kepler(y);
end
end

function y_values = bak_euler(y0, N, t_end)
%DENNA FUNKAR INTE; VAD HÄNDER?
h = t_end/N;
y = y0;
y_values = [zeros(N+1, 1),zeros(N+1, 1),zeros(N+1, 1),zeros(N+1, 1)];
for i  = 1:N
    for k = 1:4
       y_values(i,k) = y(k);
    end
    y = newton_bak(y, 10^(-6), h);
end
end

function y_values = symp_euler(y0, N, t_end)
h = t_end/N;
y = y0;
y_values = [zeros(N+1, 1),zeros(N+1, 1),zeros(N+1, 1),zeros(N+1, 1)];
for i = 1:N+1
   for k = 1:4
       y_values(i,k) = y(k);
   end
   q = [y(1);y(2)];
   p = [y(3);y(4)];
   pd = pdot(q);
   y = y + [h*p(1) + h^2*pd(1);h*p(2) + h^2*pd(2);h*pd(1);h*pd(2)];
end
end

function pd = pdot(q)
%Beräknar vektorn pdot från vektor q. 
pd = -[q(1)/(q(1)^2+q(2)^2)^(3/2);q(2)/(q(1)^2+q(2)^2)^(3/2)];
end

function y = newton_bak(y0, tol, h)
maxiter = 100;
f = @(y) y-h*kepler(y)-y0;
Df = @(y) jacobian_bak(y, h);
s = inf;
y = y0;
iter = 0;
while norm(s) >= tol
    if iter >= maxiter
        disp(iter)
        disp('Maxiterationer är nådda')
        break
    end
    s = -Df(y)\f(y);
    y = y + s;
    iter  = iter+1;
end
end

function Df = jacobian_bak(y, h)
%Beräknar jacobianen bla balb alb 
Df = zeros(4,4);
q = [y(1);y(2)];
r = norm(q);
Df(3,1) = -h/r^3*(-1+3*q(1)^2/r^2);
Df(3,2) = -h*3*q(1)*q(2)/r^5;
Df(4,1) = -h*3*q(1)*q(2)/r^5;
Df(4,2) = -h/r^3*(-1+3*q(2)^2/r^2);
Df(1,3) = -h;
Df(2,4) = -h;
Df = eye(4) + Df;
end

function H = hamiltonian(y)
%Beräknar hamiltonianen
q = [y(1);y(2)];
p = [y(3);y(4)];
H = 1/2*norm(p) - 1/(norm(q))^2;
end

function E_values = energi(y_values)
%Tar in y-värden från en lösningsmetod och beräknar energin för varje
%tidssteg
E_values = zeros(length(y_values), 1);
for i = 1:length(y_values)
    E_values(i) = hamiltonian(y_values(i, :));
end
end
