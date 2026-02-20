%% Uppgift 1
a = 0.5;
y0 = [1-a;0;0;sqrt((1+a)/(1-a))];
t_end = 100;
N1 = 200000; %Första ordningens metoder
h1 = t_end/N1;

N2 = 4000; %Andra ordningens metoder
h2 = t_end/N2;

%Framåt Euler
y_feuler = fram_euler(y0, N1, t_end);
energi_feuler = energi(y_feuler);
q1f = y_feuler(:,1);
q2f = y_feuler(:,2);

figure
plot(q1f, q2f)
title('Framåt Euler')

%Bakåt Euler
y_beuler = bak_euler(y0, N1, t_end);
energi_beuler = energi(y_beuler);
q1b = y_beuler(:,1);
q2b = y_beuler(:,2);

figure
plot(q1b, q2b)
title('Bakåt Euler')

%Implicita mittpunktsmetoden
y_mitt = imp_mitt(y0, N2, t_end);
energi_mitt = energi(y_mitt);
q1t = y_mitt(:, 1);
q2t = y_mitt(:, 2);

figure
plot(q1t, q2t)
title('Implicita mittpunktsmetoden')

%Symplektisk Euler 
y_seuler = symp_euler(y0, N1, t_end);
energi_seuler = energi(y_seuler);
q1s = y_seuler(:,1);
q2s = y_seuler(:,2);

figure
plot(q1s, q2s)
title('Symplektisk Euler')

figure 
plot(0:h1:t_end, energi_feuler, 0:h1:t_end, energi_beuler, 0:h1:t_end, energi_seuler, 0:h2:t_end, energi_mitt)
title('Energi som funktion av tiden')
legend('Framåt', 'Bakåt', 'Symplektisk', 'Trapets')

%% Funktioner 
%ATT GÖRA: Effektivare med matriser för y, snyggare med norm(q),
%startgissning för Newtons metod, varför gör bakåt Euler ett stort hopp?

function ydot = kepler(y)
%y = (q1, q2, p1, p2)
q1dot = y(3);
q2dot = y(4);
p1dot = -1*y(1)/(y(1)^2+y(2)^2)^(3/2);
p2dot = -1*y(2)/(y(1)^2+y(2)^2)^(3/2);
ydot = [q1dot;q2dot;p1dot;p2dot];
end

function y_values = fram_euler(y0, N, t_end)
h = t_end/N;
y = y0;
y_values = [zeros(N+1, 1),zeros(N+1, 1),zeros(N+1, 1),zeros(N+1, 1)]';
for i  = 1:N+1
    y_values(:, i) = y;
    y = y+h*kepler(y); %y(n+1) = y(n) +h*y(n)dot
end
y_values = y_values';
end

function y_values = bak_euler(y0, N, t_end)
%Bakåt Euler
%Denna är första ordningen så kräver mindre tiddsteg. 
h = t_end/N;
y = y0;
y_values = [zeros(N+1, 1),zeros(N+1, 1),zeros(N+1, 1),zeros(N+1, 1)]';
for i  = 1:N+1
    y_values(:, i) = y;
    y = newton_bak(y, 10^(-6), h); %Lös olinjärt ekvationssystem med Newton
end
y_values = y_values'; %Gillar att ha q1 till p2 i kolumnerna
end

function y_values = symp_euler(y0, N, t_end)
%Symplektisk Euler
h = t_end/N;
y = y0;
y_values = [zeros(N+1, 1),zeros(N+1, 1),zeros(N+1, 1),zeros(N+1, 1)]';
for i = 1:N+1
   y_values(:, i) = y;
   q = [y(1);y(2)];
   p = [y(3);y(4)];
   pd = pdot(q); %p=tidsderivata av u
   y = y + [h*p(1) + h^2*pd(1);h*p(2) + h^2*pd(2);h*pd(1);h*pd(2)];
   %Formeln för symplektisk Euler, är egentligen bara en Taylorexpansion
end
y_values = y_values';
end

function pd = pdot(q)
%Beräknar vektorn pdot från vektor q. 
pd = -[q(1)/(q(1)^2+q(2)^2)^(3/2);q(2)/(q(1)^2+q(2)^2)^(3/2)];
end

function y = newton_bak(y0, tol, h)
%Newtons metod 
maxiter = 100;
f = @(y) y-h*kepler(y)-y0; %Denna vill vi finna rötter till
Df = @(y) jacobian_bak(y, h); %Beräknar jacobianen av f
s = inf;
y = y0;
iter = 0;
while norm(s) >= tol
    if iter >= maxiter
        disp(iter)
        disp('Maxiterationer är nådda')
        break
    end
    s = -Df(y)\f(y); %Billigare att lösa system
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

function y_values = imp_mitt(y0, N, t_end)
h = t_end/N;
y = y0;
y_values = [zeros(N+1, 1),zeros(N+1, 1),zeros(N+1, 1),zeros(N+1, 1)]';
for i =1:N+1
    y_values(:, i) = y;
    y = newton_mitt(y, 10^(-6), h);
end
y_values = y_values';
end

function y = newton_mitt(y0, tol, h)
%Newtons metod för mittpunktsmetoden, enda som skiljer sig är funktionen
%och jacobianen
f = @(y) y-y0-h*kepler(0.5*(y0+y));
Df = @(y) jacobian_mitt(0.5*(y0+y), h);
s = inf;
y = y0;
while norm(s) >= tol
    s = -Df(y)\f(y);
    y = y + s;
end
end

function Df = jacobian_mitt(y, h)
Df = zeros(4,4);
q = [y(1);y(2)];
r = norm(q);
%Väldigt lik Jacobianen för Newton Bak, 
Df(3,1) = -h/r^3*(-1+3*q(1)^2/r^2);
Df(3,2) = -h*3*q(1)*q(2)/r^5;
Df(4,1) = -h*3*q(1)*q(2)/r^5;
Df(4,2) = -h/r^3*(-1+3*q(2)^2/r^2);
Df(1,3) = -h;
Df(2,4) = -h;
Df = eye(4) + 0.5*Df; %Enda som skiljer är faktor 0.5, kedjeregeln
end

function H = hamiltonian(y)
%Beräknar hamiltonianen
q = [y(1);y(2)];
p = [y(3);y(4)];
H = 1/2*norm(p)^2 - 1/(norm(q));
end

function E_values = energi(y_values)
%Tar in y-värden från en lösningsmetod och beräknar energin för varje
%tidssteg
E_values = zeros(length(y_values), 1);
for i = 1:length(y_values)
    E_values(i) = hamiltonian(y_values(i, :));
end
end