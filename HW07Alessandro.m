clear;
u = 0.01;                   % step unit;
L = 0:u:12;                 % ft;
v0 = 300;                   % lbf;
m0 = 0;                     % lbf*ft;
w = 12.5 .* L;              % lbf/ft^2; for analytical method;
fun = @(t) 12.5 * t;        % lbf/ft^2; for numerical method;

%% Analytically determine Shear Force (Figure 2)

va = zeros(0);              % initializes empty array for shear forces V;
s = 1;                      % first position of array;
syms x                      % initializes symbolic variable;
for i = L
a = int(12.5 * x,x,0,i);    % lbf. calculates the definite integral of the 
                            % distributed load for the beam with
                            % respect to its symbolic variable x, from
                            % 0 to i;
b = double(a);              % turns syms into double;
va(s) = v0 - b;             % creates array of values of shear force;
s = s + 1;                  % moves to next position in the vector;
end

%% Numerically determine Shear Force (Figure 2)

vn= zeros(0);               % initializes empty array for shear forces V;
s = 1;                      % first position of array;

for i = L
c = integral(fun,0,i);      % lbf. Numerically approximates the integral
                            % of function "fun" from 0 to i;
vn(s) = v0 - c;             % creates array of values of shear force;
s = s + 1;                  % moves to next position in the vector;
end

%% Numerically determine Bending momentum (Figure 3)
intv = 0;                   % integral of the shear force V at 0 = 0;
mn = zeros(0);              % initializes empty array for bending momentums M;
for i = 1: (length(vn)-1)
    intv = intv + u * (( va(i)+vn(i+1) )/2);% calculates integral 
                                            % numerically as the area under
                                            % a trapezoid;
    mn(i) = m0 + intv;                      % creates array of values of
                                            % bending momentum;
end

%% Analytically determine Bending momentum (Figure 3)

ma = zeros(0);                  % initializes empty array for shear forces V;
s = 1;                          % first position of array;
syms x;                         % initializes symbolic variable;
a1 = int(12.5 * x);             % lbf; calculates indefinite integral of
                                % distributed load
for i = L
e = int( (v0 - a1) ,x,0,i );    % lbf*ft; calculates the definite integral
                                % of the shear force (a1, the indefinite 
                                % integral of distributed load for the
                                % beam) with respect to its symbolic
                                % variable x, from 0 to i;
f = double(e);                  % turns syms into double;
ma(s) = m0 + f;                 % creates array of values of shear force;
s = s + 1;                      % moves to next position in the vector;
end;

%% Plot Internal bending moment (3)
figure(3)                               
set(gcf, 'units','normalized','outerposition',[0.66, 0.3,  0.33, 0.5]);
plot(L,ma,'g','LineWidth',10);              % Analytical
hold on
plot(L(1:(end-1)),mn,'b','LineWidth',5);    % Numerical
xlabel('Position (ft)');
ylabel('Internal bending moment (lbf*ft)');
legend('Analytical','Numerical','Location','SE');
% area trapezoid = ( (a + b)/2 ) * h;
%% Plot Internal shear force (2)
figure(2)
set(gcf,'units','normalized','outerposition',[0.33, 0.3, 0.33, 0.5]);
plot(L,va,'g','LineWidth',10)               % Analytical
xlabel('Position (ft)');
ylabel('Internal shear force (lbf)');
hold on

plot(L,vn,'b','LineWidth',5)                % Numerical
xlabel('Position (ft)');
ylabel('Internal shear force (lbf)');
hold on

D = diff(mn)./diff(L(1:(end-1)));           % differentiates numerically
plot(L(1:(end-2)),D,'y','LineWidth',1);     % Reconstructed
legend('Analytical','Numerical','Reconstructed','Location','NE');

%% Plot Distributed load (1)
figure(1);
set(gcf,'units','normalized','outerposition',[0.001, 0.3, 0.33, 0.5]);
plot(L,w,'b','LineWidth',5);                % Numerical
xlabel('Position (ft)');
ylabel('Distributed load (lbf/ft)');
hold on

DD = abs( diff(D)./diff(L(1:(end-2))) );    % Differentiates numerically
plot(L(1:(end-3)),DD,'y','LineWidth',1);    % Reconstructed
legend('Original','Reconstructed','Location','SE')
