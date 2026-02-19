%% PARAMETRAR
N=200; % antal intervall
T=2; % sluttid
dx=1/N; % steglängd i rummet
dt=dx/2.0; % tidssteg, tänk på stabilitetsvillkoren
M=round(T/dt); % antal tidsteg
c=1; % våghastighet
% allokering av minne
u=zeros(N-1,M+1); % u(n,m) lösningens värde vid tid (m-1)*dt i position n*dx
p=zeros(N-1,M+1); % p=u’
A=zeros(N-1,N-1); % Au är differensapproximation av d^2 u/dx^2
x = dx*(1:N-1)'; % x(n) är n*dx
E = zeros(1,M+1); % För att beräkna energin i varje tidssteg.

%% Skapa matrisen A
A = A -2*eye(N-1);
B = zeros(N-1,N-1);
for i = 1:N-2
    B(i, i+1) = 1;
end
A = 1/(dx)^2*(A +B +B');
%% Sätt begynnelsedata för u och p.
g  = @(x) exp(-200*(x-0.5)^2);
for i = 1:N-1
    u(i,1) = g(i*dx);
end
%Händer inget med p
%% Räkna ut energin E vid tiden 0.
E(1, 1) = energy(u(:, 1), p(:, 1), c, A);
nframe=M+1; % kommando för film
mov(1:nframe)= struct('cdata',[],'colormap',[]);
figure;
plot(x,u(:,1), 'b', 'Linewidth', 1); %Plot vid tiden t=0.
axis([0 1 -1 1])
set(gca, 'nextplot', 'replacechildren')
drawnow
mov(1)=getframe(gcf); %Första frame i filmen.

%% tidstegning med symplektiska Euler
L=1;
for m=1:M 
p(:, m+1) = p(:, m) + dt*A*u(:,m);
u(:, m+1) = u(:, m) + dt*p(:,m+1) + dt^2*A*u(:,m);

X = [0;x;L]; U = [0;u(:,m);0];
plot(X, U, 'b', 'Linewidth', 1)
hold on;
t = m*dt;

%%Plotta även lösningen från d’Alemberts formel
ua = u;
for j = 1:N-1
    ua(j, m+1) = 0.5*(g(x(j)+c*t)-g(x(j)-c*t));
end

%Xa = [0;x;L]; Ua = [0;ua(:,m);0];
%plot(Xa, Ua, X, U, 'b', 'Linewidth', 1)
%hold on;

text(0.05,-0.8, sprintf('t=%.2f', t))
set(gca, 'nextplot', 'replacechildren')
drawnow
mov(m+1)=getframe(gcf);

%Räkna ut energin av den numeriska lösningen vid detta tidsstag
E(m+1) = energy(u(:, m+1), p(:, m+1), c, A);
end

function en  = energy(u, p, c, A)
%tar in vektorer u och p och beräknar energin för en tidpunkt
udd = c^2*A*u;
en = 0.5*(norm(p)^2 -(dot(u, udd)));
end