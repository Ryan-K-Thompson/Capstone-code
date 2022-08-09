function MakeNozzle(AR)
% Make nozzle completely user-defined.
% Variables
% Length
% Nozzle area ratio
% if nargin == 0
%     NL = input('Enter the nozzle length: ');
%     AR = input('Enter the nozzle area ratio: ');
% else
%     NL = varargin{1,1};
%     AR = varargin{1,2};
% end
NL = 10;
% Initialize coordinates
x = (0:1e-3:NL)';
L = length(x);
% Reservoir should be at 1/6 of total nozzle length
% Throat should be at 2/5 of total nozzle length
% Nozzle shape is flat from start to xR.
% Follow a downslope from the next element to xR
tol = 1e-3;
for i = 1:L
    if abs(x(i)-(NL/8)) < tol
        xRct = i;
    elseif abs(x(i)-(0.5*NL)) < tol
        xTHct = i;
        break
    end
end
NW = 10;
Constant = nan(L,1);
Converging = Constant;
Diverging = Converging;
TA = NW/3;
Diverging(xTHct) = TA;
Diverging(xRct+1) = NW;
% Construct converging cosine wave
ampC = (NW-TA)/2;
wC = (pi)/(x(xTHct)-x(xRct));
psC = -wC*x(xRct);
CW = @(x) ampC*cos(wC*x+psC)+TA+ampC;
% Construct diverging cosine wave
EA = AR*TA;
ampD = (EA-TA)/2;
wD = pi/(x(end)-x(xTHct));
psD = -wD*x(end);
DW = @(x) ampD*cos(wD*x+psD)+TA+ampD;
% Make nozzle
Constant(1:xRct,1) = NW*ones(xRct,1);
for i = 1:L
    if i > xRct && i < xTHct
        Converging(i) = CW(x(i));
    elseif i > xTHct
        Diverging(i) = DW(x(i));
    end
end
figure
plot(x,Constant,x,Converging,x,-Constant,x,-Converging, ... 
    x,Diverging,x,-Diverging,'Color','b','LineWidth',4);
grid on
if max(Constant) > max(Diverging)
    M = max(Constant);
else
    M = max(Diverging);
end
T = sprintf('User-Defined Nozzle AR = %4.2f',AR);
title(T);
axis([x(1) 1.3*x(end) 1.2*-M 1.2*M]);

end