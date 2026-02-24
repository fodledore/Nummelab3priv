%uppgift 3

N2 = 20000;

a = 0.5;
y0 = [1-a;0;0;sqrt((1+a)/(1-a))];
t_end = 100;
tspan = linspace(0, t_end, N2);

% ode45
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
tic
[t_ode, y_ode] = ode45(@kepler, tspan, y0, options);
time_ode = toc;
energi_ode = energy(y_ode);
q1ode = y_ode(:,1);
q2ode = y_ode(:,2);
y_ref = y_ode(end,:); % y-vektorn vid t = 100, används som referens

figure
plot(q1ode, q2ode)
title('ode45')

% symplektisk euler
y_seuler = symp_euler(y0, N2, t_end);
energi_seuler = energy(y_seuler);
q1s = y_seuler(:,1);
q2s = y_seuler(:,2);

figure
plot(q1s, q2s)
title('Symplektisk Euler')

% Implicita mittpunktsmetoden
y_mitt = imp_mitt(y0, N2, t_end);
energi_mitt = energy(y_mitt);
q1m = y_mitt(:, 1);
q2m = y_mitt(:, 2);

figure
plot(q1m, q2m)
title('Implicita mittpunktsmetoden')


% ---------------------------- (fel i lösningsvektorerna)


N_vals = [1300 1700 3000 4000 5000 6000 8000 10000 12000 14000 16000 20000]; % vi testar olika N-värden
errors_se = zeros(size(N_vals));
errors_mitt = zeros(size(N_vals));
times_se = zeros(size(N_vals));
times_mitt = zeros(size(N_vals));
h_vals = t_end ./ N_vals;

for k = 1:length(N_vals)
    N = N_vals(k);
    
    % Symplektisk Euler
    tic
    y_se = symp_euler(y0, N, t_end);
    times_se(k) = toc;
    errors_se(k) = norm(y_se(end,:) - y_ref);
    
    % Implicit midpoint
    tic
    y_mittp = imp_mitt(y0, N, t_end);
    times_mitt(k) = toc;
    errors_mitt(k) = norm(y_mittp(end,:) - y_ref);
end

% --------------------- (plotta fel mot steglängd)

figure
loglog(h_vals, errors_se, 'o-', 'DisplayName','Symplektisk Euler')
hold on
loglog(h_vals, errors_mitt, 's-', 'DisplayName','Implicit midpoint')
xlabel('h')
ylabel('Fel i lösningsvektor vid t = 100')
title('Fel vs steglängd')
legend
grid on

% --------------- (beräkna uppmätt noggranhetsordning)

p_se = log(errors_se(1:end-1)./errors_se(2:end)) ./ log(h_vals(1:end-1)./h_vals(2:end));
p_mitt = log(errors_mitt(1:end-1)./errors_mitt(2:end)) ./ log(h_vals(1:end-1)./h_vals(2:end));

disp(['Uppmätt ordning symplektisk Euler: ', num2str(mean(p_se))])
disp(['Uppmätt ordning implicit midpoint: ', num2str(mean(p_mitt))])

% ------------------ (plotta beräkningstid mot fel)

figure
loglog(errors_se, times_se, 'o-', 'DisplayName','Symplektisk Euler')
hold on
loglog(errors_mitt, times_mitt, 's-', 'DisplayName','Implicit midpoint')
xlabel('Fel i lösningsvektor')
ylabel('Beräkningstid [s]')
title('Effektivitet: tid vs fel')
legend
grid on

% ------------------- (plotta fel i energi över tid)

y_seuler = symp_euler(y0, N_vals(end), t_end);
y_mitt = imp_mitt(y0, N_vals(end), t_end);

E0 = energy(y0);
E_se = zeros(N_vals(end)+1,1);
E_mitt = zeros(N_vals(end)+1,1);
%E_ode = zeros(length(t_ode),1);

for i = 1:N_vals(end)+1
    E_se(i) = energy(y_seuler(i,:)');
    E_mitt(i) = energy(y_mitt(i,:)');
    %E_ode(i) = energy(y_ode(i,:)');
end

figure
plot(abs(E_se - E0), 'b', 'DisplayName', 'Symplektisk Euler')
hold on
plot(abs(E_mitt - E0), 'r', 'DisplayName', 'Implicit midpoint')
%hold on
%plot(abs(E_ode - E0), 'r', 'DisplayName', 'ODE45')
xlabel('Steg')
ylabel('|E - E_0|')
title('Energifel över tid: jämförelse av metoder')
legend
grid on

% ----------- (plottar maxfel i energi för olika N-värden)

maxE_se = zeros(size(N_vals));
maxE_mitt = zeros(size(N_vals));

E0 = energy(y0);   % exakt energi (konstant)

for k = 1:length(N_vals)
    
    N = N_vals(k);
   
    y_se = symp_euler(y0, N, t_end);
    E_se = zeros(N+1,1);
    
    for i = 1:N+1
        E_se(i) = energy(y_se(i,:)');
    end
    
    maxE_se(k) = max(abs(E_se - E0));
    
    y_mitt = imp_mitt(y0, N, t_end);
    E_mitt = zeros(N+1,1);
    
    for i = 1:N+1
        E_mitt(i) = energy(y_mitt(i,:)');
    end
    
    maxE_mitt(k) = max(abs(E_mitt - E0));
    
end

figure
loglog(h_vals, maxE_se, 'o-', 'DisplayName','Symplektisk Euler')
hold on
loglog(h_vals, maxE_mitt, 's-', 'DisplayName','Implicit midpoint')
xlabel('h')
ylabel('Maximalt energifel')
title('Max energifel vs steglängd')
legend
grid on


% ----------------- (fel i lösningsvektor vid sista tidpunkten, t=100)

err_se = norm(y_seuler(end,:) - y_ref);
err_mitt = norm(y_mitt(end,:) - y_ref);
disp(['Fel symplektisk Euler (N=', num2str(N_vals(end)), '): ', num2str(err_se)])
disp(['Fel implicit midpoint (N=', num2str(N_vals(end)), '): ', num2str(err_mitt)])

% ------------------- (ODE45 – separat analys)

tolerances = [1e-4 1e-6 1e-8 1e-10 1e-12];

errors_ode = zeros(size(tolerances));
times_ode = zeros(size(tolerances));
maxE_ode = zeros(size(tolerances));

E0 = energy(y0);

for k = 1:length(tolerances)

    options = odeset('RelTol',tolerances(k),'AbsTol',tolerances(k));

    tic
    [t_tmp, y_tmp] = ode45(@kepler, [0 t_end], y0, options);
    times_ode(k) = toc;

    % fel i lösningsvektor vid t = 100
    errors_ode(k) = norm(y_tmp(end,:) - y_ref);

    % energifel
    E_tmp = zeros(length(t_tmp),1);
    for i = 1:length(t_tmp)
        E_tmp(i) = energy(y_tmp(i,:)');
    end

    maxE_ode(k) = max(abs(E_tmp - E0));

    % spara en körning (den strängaste toleransen) för energifel vs tid
    if k == length(tolerances)
        t_energy = t_tmp;
        E_energy = abs(E_tmp - E0);
    end

end


% -------------------(fel mot tolerans för ODE45)

figure
loglog(tolerances, errors_ode, 'o-')
xlabel('Tolerans')
ylabel('Fel i lösningsvektor vid t = 100')
title('ode45: Fel vs tolerans')
grid on

% ---------------- (Effektivitet: beräkningstid mot fel för ODE45)

figure
loglog(errors_ode, times_ode, 'o-')
xlabel('Fel i lösningsvektor')
ylabel('Beräkningstid [s]')
title('ode45: Effektivitet (tid vs fel)')
grid on

% --------------------(Energifel mot tid för ODE45) 

figure
plot(t_energy, E_energy)
xlabel('Tid')
ylabel('|E - E_0|')
title('ode45: Energifel över tid')
grid on


% ---------------(Max energifel mot tolerans för ODE45)

figure
loglog(tolerances, maxE_ode, 'o-')
xlabel('Tolerans')
ylabel('Maximalt energifel')
title('ode45: Max energifel vs tolerans')
grid on

% =========================(funktioner)


% ger derivatan till y
function ydot = kepler(~, y)
    %y = (q1, q2, p1, p2)
    q1dot = y(3);
    q2dot = y(4);
    p1dot = -1*y(1)/(y(1)^2+y(2)^2)^(3/2);
    p2dot = -1*y(2)/(y(1)^2+y(2)^2)^(3/2);
    ydot = [q1dot;q2dot;p1dot;p2dot];
end


% symplektisk euler
function y_values = symp_euler(y0, N, t_end)
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

% Beräknar vektorn pdot från vektor q
function pd = pdot(q)
    pd = -[q(1)/(q(1)^2+q(2)^2)^(3/2);q(2)/(q(1)^2+q(2)^2)^(3/2)];
end


% beräknar energi
function E = energy(y)
    q1 = y(1); q2 = y(2);
    p1 = y(3); p2 = y(4);
    E = 0.5*(p1^2+p2^2) - 1/sqrt(q1^2+q2^2);
end


% implicita mittpunktsmetoden
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
    f = @(y) y-y0-h*kepler(0, 0.5*(y0+y));
    Df = @(y) jacobian_trap(0.5*(y0+y), h);
    s = inf;
    y = y0;
    while norm(s) >= tol
        s = -Df(y)\f(y);
        y = y + s;
    end
end

function Df = jacobian_trap(y, h)
    Df = zeros(4,4);
    q = [y(1);y(2)];
    r = norm(q);
    %Väldigt lik Jacobianen för Newton Bak, endast skalad med 0.5
    Df(3,1) = -h/r^3*(-1+3*q(1)^2/r^2);
    Df(3,2) = -h*3*q(1)*q(2)/r^5;
    Df(4,1) = -h*3*q(1)*q(2)/r^5;
    Df(4,2) = -h/r^3*(-1+3*q(2)^2/r^2);
    Df(1,3) = -h;
    Df(2,4) = -h;
    Df =(eye(4) + 0.5*Df);
end