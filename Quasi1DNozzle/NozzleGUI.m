function NozzleGUI
%% Information
% Kazuaki Iida
% Updated 9/1/14

%% Application Code
close all force
global g NPR AR AreaRatio rmin radius

% Construct the components
GUI = figure('Visible','off','Position',[200,100,1000,600]);
hlabel = uicontrol('Style','text','String','Nozzle GUI', ... 
    'Position',[330,560,250,35]);
set(hlabel,'FontSize',22);
instruct = sprintf('Enter the nozzle pressure ratio, gamma, and area ratio.');
hlabel2 = uicontrol('Style','text','String',instruct, ... 
    'Position',[710,445,250,88]);
set(hlabel2,'FontSize',18);
hlabel2nd = uicontrol('Style','text','String','Nozzle Pressure Ratio','Position',[700,390,180,20]);
set(hlabel2nd,'FontSize',12);
hedit = uicontrol('Style','edit','Position',[910,390,40,20],'Callback',@Activate);
set(hedit,'FontSize',12);
hlabel3rd = uicontrol('Style','text','String','Gamma','Position',[730,330,100,20]);
set(hlabel3rd,'FontSize',12);
hlabelAR = uicontrol('Style','text','String','Area Ratio','Position',[730,360,100,20]);
set(hlabelAR,'FontSize',12);
heditAR = uicontrol('Style','edit','Position',[910,360,40,20],'Callback',@Activate3);
set(heditAR,'FontSize',12);
hedit2 = uicontrol('Style','edit','Position',[910,330,40,20],'Callback',@Activate2);
set(hedit2,'FontSize',12);
push2 = uicontrol('Style','pushbutton','String','Plot','Position',[850,220,100,60],'Callback',@Plot);
set(push2,'FontSize',14);
pushRestart = uicontrol('Style','pushbutton','String','Restart GUI', ... 
    'Position',[850,140,120,60],'Callback',@Restart);
set(pushRestart,'FontSize',14);
push = uicontrol('Style','pushbutton','String','Exit GUI', ... 
    'Position',[850,60,100,60],'Callback',@Button);
set(push,'FontSize',14);

% Add spaces for plots
hpNozzle = axes('Units','Pixels','Position',[340,320,330,210]);
hpMach = axes('Units','Pixels','Position',[60,60,330,210]);
hpPressure = axes('Units','Pixels','Position',[460,60,330,210]);

% Make GUI visible
set(GUI,'Visible','on');
set(gcf,'Color',[0.68,0.68,0.68]); % Makes figure gray

% Initialize NPR and g to nan
NPR = nan;
g = nan;
AR = nan;

% Intialize Counter
P = 0;

% Get NPR
    function Activate(source1,eventdata,handles)
        NPR = get(source1,'String');
        NPR = str2double(NPR);
    end

% Get gamma
    function Activate2(source2,eventdata,handles)
        g = get(source2,'String');
        g = str2double(g);
    end

    function Activate3(source,eventdata,handles)
        AR = str2double(get(source,'String'));
    end

    function Button(source,eventdata)
        fprintf('Closing GUI.\n');
        close all force
    end

% Plot
    function Plot(hObject,eventdata,handles,varargin)
        if NPR < 1
           fprintf('\nNozzle pressure ratio must be greater than 1.\n');
           NozzleGUI;
        elseif g == 0 && NPR == 0
            fprintf('\nGamma must be greater than zero.\n');
            NozzleGUI;
            return
        elseif isnan(g) == 1 && isnan(NPR) == 1 && isnan(AR) == 1
            clc
            fprintf('\nEnter real values.\n');
            NozzleGUI;
            return
        elseif AR < 1
            % AR is defined as the exit area/throat area, so if AR is less 
            % than 1, the exit would be the throat.
            fprintf('Area ratio must be greater than or equal to 1.\n');
            fprintf('Calculation made based on AR = 1.\n');
            AR = 1;
            quasi1dflow(NPR,g,AR);
        else
            quasi1dflow(NPR,g,AR);
            P = P+1;
        end
    end

    function Restart(hObject,eventdata)
        clc
        NozzleGUI;
    end


function quasi1dflow(NPR,g,AR)
% Map 

% Ask for nozzle pressure ratio 

% Depending on the value of the nozzle pressure ratio, 
% there are 5 cases. 

% Case 1: Flow is not choked. Subsonic flow everywhere. 
% Depending on how close the pressure ratio is to the subsonic solution, 
% the Mach number at the throat should be closer to 1. 

% Case 2: Subsonic and choked. Isentropic. 

% Case 3: NPR is between subsonic and supersonic solution. 
% Either a normal shock exists inside the nozzle, or an oblique shock 
% exists outside the nozzle. 
% If a normal shock is found inside the nozzle, draw the location of the
% shock. 
% If an oblique shock exists outside the nozzle, draw the orienation of
% the shock depending on how close the pressure ratio is to the supersonic 
% solution. 

% Case 4: Perfectly expanded. 

% Case 5: Underexpanded jet. Draw Expansion fans outside the nozzle. 
% Pressure and Mach number outside the nozzle change according ot Prandtl
% Meyer expansion formula

% Defined NPR as the inverse, so fix it. 
NPR = 1/NPR;

AreaRatio = AR;

NL = 2;
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
Reservoir= nan(L,1);
Converging = Reservoir;
Diverging = Reservoir;
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
Reservoir(1:xRct,1) = NW*ones(xRct,1);
for i = 1:L
    if i > xRct && i < xTHct
        Converging(i) = CW(x(i));
    elseif i > xTHct
        Diverging(i) = DW(x(i));
    end
end
radius = zeros(L,1);
for i = 1:L
    if isnan(Reservoir(i)) == 0
        radius(i) = Reservoir(i);
    elseif isnan(Converging(i)) == 0
        radius(i) = Converging(i);
    elseif isnan(Diverging(i)) == 0
        radius(i) = Diverging(i);
    end
end
rmin = min(radius);
for i = 1:L
    if radius(i) == rmin
        Count = i;
        break
    end
end
Mach = nan(L,1);
M = Mach;
for i = 1:L
    f = @(M) (1/M^2)*((2/(g+1))*(1+((g-1)/2)*M^2))^((g+1)/(g-1))... 
       - (radius(i)/rmin)^2;
    if i < Count
       Mach(i) = bisect(f,0.01,1,100);
       M(i) = Mach(i);
    elseif i == Count
       Mach(i) = 1;
       M(i) = 1;
    else
       Mach(i) = bisect(f,1,10,100);
       M(i) = bisect(f,0.01,1,100);
    end       
end
PresRatio = nan(L,1);
% Pressure Ratio for subsonic choked solution
for i = 1:L
    PresRatio(i) = PR(M(i));
end
PresRatio2 = nan(L,1);
% Pressure Ratio for supersonic solution (perfectly expanded)
for i = 1:L
    PresRatio2(i) = PR(Mach(i));
end
Subsonic = PresRatio(end);
Supersonic = PresRatio2(end);    

    function PlotNozzle(varargin)
        Parameter = varargin{1,1};
        if Parameter == 1
            % Normal shock
            SL = varargin{1,2};
            index = varargin{1,3};
        hNozzle = plot(hpNozzle,x,Reservoir, ... 
            x,Converging, ... 
            x,-Reservoir, ... 
            x,-Converging, ... 
            x,Diverging, ... 
            x,-Diverging, ...
            x(Count)*ones(1,10),linspace(-rmin,rmin,10),SL*ones(1,10),linspace(-radius(index),radius(index),10));
        for NozzleIndex = 1:length(hNozzle)-2
                set(hNozzle(NozzleIndex),'LineWidth',3,'Color','b');
        end
            set(hNozzle(length(hNozzle)-1),'LineWidth',1.2,'Color','k');
            set(hNozzle(length(hNozzle)),'LineWidth',5,'Color','r');
        elseif Parameter == 2
            % Expansion fan
            MachEnd=varargin{1,3};
            MachOut=varargin{1,4};
            EFx = (x(end):0.01:2.5)';
            EF = nan(length(EFx),1);
            EF2 = EF;
            EF3 = EF;
            % slope of the tail is the angle
            angle = varargin{1,2};
            slope1 = -tan(asin(1/MachEnd));
            slope3 = tan(-asin(1/MachEnd)+angle);
            slope2=0;
            EF(1)=Diverging(end);
            EF2(1)=EF(1);
            EF3(1)=EF(1);
            for j = 2:length(EFx)
                EF(j)=EF(j-1)+slope1*(EFx(j)-EFx(j-1));
                EF2(j)=EF2(j-1)+slope2*(EFx(j)-EFx(j-1));
                EF3(j)=EF3(j-1)+slope3*(EFx(j)-EFx(j-1));
            end
        hNozzle = plot(hpNozzle,x,Reservoir, ... 
            x,Converging, ... 
            x,-Reservoir, ... 
            x,-Converging, ... 
            x,Diverging, ... 
            x,-Diverging, ...
            x(Count)*ones(1,10),linspace(-rmin,rmin,10), ...
            EFx,EF, ...
            EFx,-EF, ...
            EFx,EF2, ...
            EFx,-EF2, ...
            EFx,EF3, ...
            EFx,-EF3);           
        for NozzleIndex = 1:length(hNozzle)
            if NozzleIndex <= 6
                set(hNozzle(NozzleIndex),'LineWidth',3,'Color','b');
            elseif NozzleIndex == 7
                set(hNozzle(NozzleIndex),'LineWidth',1.2,'Color','k');
            else
                set(hNozzle(NozzleIndex),'LineWidth',2,'Color','g')
            end
        end
        elseif Parameter == 3
            % Oblique shock
            angle = varargin{1,2};
            OSx = (x(end):0.01:2.4)';
            OS = nan(length(OSx),1);   
            OS(1) = Diverging(end);
            slope = tan(angle);
            for OSIndex = 2:length(OSx)
                OS(OSIndex) = OS(OSIndex-1)-slope*(OSx(OSIndex)-OSx(OSIndex-1));
            end
            hNozzle = plot(hpNozzle,x,Reservoir, ... 
            x,Converging, ... 
            x,-Reservoir, ... 
            x,-Converging, ... 
            x,Diverging, ... 
            x,-Diverging, ...
            x(Count)*ones(1,10),linspace(-rmin,rmin,10), ...
            OSx,OS,OSx,-OS);
            for NozzleIndex = 1:length(hNozzle)
                if NozzleIndex <= 6
                    set(hNozzle(NozzleIndex),'LineWidth',3,'Color','b');
                elseif NozzleIndex == 7
                    set(hNozzle(NozzleIndex),'LineWidth',1.2,'Color','k');
                else
                    set(hNozzle(NozzleIndex),'LineWidth',5,'Color','r');
                end
            end
        else
            % Nothing (no shocks or expansion fans)
        hNozzle = plot(hpNozzle,x,Reservoir, ... 
            x,Converging, ... 
            x,-Reservoir, ... 
            x,-Converging, ... 
            x,Diverging, ... 
            x,-Diverging, ...
            x(Count)*ones(1,10),linspace(-rmin,rmin,10));      
        for NozzleIndex = 1:length(hNozzle)-1
            set(hNozzle(NozzleIndex),'LineWidth',3,'Color','b');
        end
            set(hNozzle(length(hNozzle)),'LineWidth',1.2,'Color','k');        
        end
        T = sprintf('CONVERGING-DIVERGING NOZZLE for NPR = %4.3f',1/NPR);
        title(hpNozzle,T);
        if radius(end) > 10
            set(hpNozzle,'XLim',[x(1) 1.5*x(end)],'YLim',[-1.5*radius(end) 1.5*radius(end)]);
        else
            set(hpNozzle,'XLim',[x(1) 1.5*x(end)],'YLim',[-11 11]);
        end
        hThroat = uicontrol('Style','text','String','Thin black line is the throat','Position',[20,500,200,20]);
        set(hThroat,'FontSize',12);
        hShock  = uicontrol('Style','text','String','Red lines are shocks','Position',[20,450,200,20]);
        set(hShock,'FontSize',12);
        hEF = uicontrol('Style','text','String','Green lines are expansion fans','Position',[20,390,200,40]);
        set(hEF,'FontSize',12);
    end

    function varargout = PlotMach(Parameter,NPR,Case)
        % Only use input if the solution does not correspond to 
        % either the subsonic or supersonic solution. 
        % Parameter contains the Mach number to be plotted 
        % (if the solution contains shocks or expansion fans, or 
        % is subsonic everywhere), OR it contains the scalar value
        % 0 or 1, indicating the subsonic or supersonic solution.
        % Case = 1: Shock
        % Case = 2: Expansion fan
        % Case = 0: Neither
        if length(Parameter) ~= 1 && Case == 1
            % Shock
            % Feed in shock location and critical Mach number (pre-shock
            % Mach number)
            xnew = [x;(x(end)+0.001:0.001:3)'];
            Z = length(xnew);
            Parameter = [Parameter;Parameter(end)*ones(Z-L,1)];
            hMach = plot(hpMach,xnew,Parameter,x(end)*ones(10,1),(linspace(min(Parameter),1.1*max(Parameter),10))');
            set(hMach(1),'LineWidth',3,'Color','b');
            set(hMach(2),'LineWidth',2,'Color','k');
        elseif length(Parameter) ~= 1 && Case == 2
            % Oblique shock
            % Feed in supersonic solution and exit Mach number
            MachExit = NPR;
            xnew = [x;(x(end)+0.001:0.001:3)'];
            Z = length(xnew);
            Parameter = [Parameter;Parameter(end)*ones(100,1);MachExit*ones(Z-L-100,1)];
            hMach = plot(hpMach,xnew,Parameter,x(end)*ones(10,1),(linspace(min(Parameter),1.1*max(Parameter),10))');
            set(hMach(1),'LineWidth',3,'Color','b');
            set(hMach(2),'LineWidth',2,'Color','k');
        elseif length(Parameter) == 1 && Case == 2
            % Expansion fan
            % Expansion fans are gradual, so the solution should be 
            % a smooth curve. 
            % Use a cubic spline interpolant.
            F = @(M2) ((1+((g-1)/2)*Mach(end)^2)/(1+((g-1)/2)*M2^2))^(g/(g-1)) - (NPR/PresRatio2(end));
            M2 = fzero(F,1.1*Mach(end));
            if M2 < 0
                M2 = -M2;
            end
            varargout{1,1} = M2;
            xnew = [x;(x(end)+0.001:0.001:3)'];
            Z = length(xnew);
            LSpline = 400;
            S1 = [x(end-10:end);xnew(end-LSpline:end)];
            S2 = [Mach(end-10:end);M2*ones(LSpline+1,1)];
            A = nan(Z-L,1);
            for z = 1:length(A)
                A(z) = spline(S1,S2,xnew(L+z));
            end
            Machnew = [Mach;A];
            hMach = plot(hpMach,xnew,Machnew,x(end)*ones(10,1),(linspace(min(Machnew),1.1*max(Machnew),10))');
            set(hMach(1),'LineWidth',3,'Color','b');
            set(hMach(2),'LineWidth',2,'Color','k');
        elseif length(Parameter) ~= 1 && Case == 0
            % Subsonic but not choked
            xnew = [x;(x(end)+0.001:0.001:3)'];
            Z = length(xnew);
            Parameter = [Parameter;Parameter(end)*ones(Z-L,1)];
            hMach = plot(hpMach,xnew,Parameter,x(end)*ones(10,1),(linspace(min(Parameter),1.1*max(Parameter),10))');
            set(hMach(1),'LineWidth',3,'Color','b');
            set(hMach(2),'LineWidth',2,'Color','k');
        elseif Parameter == true
            % Plot supersonic solution
            xnew = [x;(x(end)+0.001:0.001:3)'];
            Z = length(xnew);
            Machnew = [Mach;Mach(end)*ones(Z-L,1)];
            hMach = plot(hpMach,xnew,Machnew,x(end)*ones(10,1),(linspace(min(Mach),1.1*max(Mach),10))');
            set(hMach(1),'LineWidth',3,'Color','b');
            set(hMach(2),'LineWidth',2,'Color','k');
        elseif Case == -1
            % No flow
            xnew = [x;(x(end)+0.001:0.001:3)'];
            Z = length(xnew);
            MachZero = zeros(Z,1);
            hMach = plot(hpMach,xnew,MachZero,x(end)*ones(10,1),(linspace(0,0.01,10))');
            set(hMach(1),'LineWidth',3,'Color','b');
            set(hMach(2),'LineWidth',2,'Color','k');
        else
            % Plot subsonic solution 
            xnew = [x;(x(end)+0.001:0.001:3)'];
            Z = length(xnew);
            Mnew = [M;M(end)*ones(Z-L,1)];
            hMach = plot(hpMach,xnew,Mnew,x(end)*ones(10,1),(linspace(min(M),1.1*max(M),10))');
            set(hMach(1),'LineWidth',3,'Color','b');
            set(hMach(2),'LineWidth',2,'Color','k');      
        end
        title(hpMach,'Mach number plot');
        legend(hpMach,'Current Solution','Nozzle exit','Location','Best');
        set(get(hpMach,'XLabel'),'String','Horizontal Coordinate');
        set(get(hpMach,'YLabel'),'String','Mach number');
        grid(hpMach);
    end

    function PlotPressure(Parameter,NPR,Case)
        % Only use input if the solution does not correspond to 
        % either the subsonic or supersonic solution. 
        % Parameter contains the Mach number to be plotted 
        % (if the solution contains shocks or expansion fans, or 
        % is subsonic everywhere), OR it contains the scalar value
        % 0 or 1, indicating the subsonic or supersonic solution.
        if length(Parameter) ~= 1 && Case == 1
            % Shock
            % Feed in pressure ratio
            xnew = [x;(x(end)+0.001:0.001:3)'];
            Z = length(xnew);
            Parameter = [Parameter;Parameter(end)*ones(Z-L,1)];
            hPressure = plot(hpPressure,xnew,Parameter,x(end)*ones(10,1),(linspace(0,1.1,10))');
            set(hPressure(1),'LineWidth',3,'Color','b');
            set(hPressure(2),'LineWidth',2,'Color','k');
        elseif length(Parameter) ~= 1 && Case == 2
            % Oblique shock
            xnew = [x;(x(end)+0.001:0.001:3)'];
            Z = length(xnew);
            % Feed in supersonic solution 
            % Have shock away from the nozzle exit
            Parameter = [Parameter;Parameter(end)*ones(100,1);NPR*ones(Z-L-100,1)];
            hPressure = plot(hpPressure,xnew,Parameter,x(end)*ones(10,1),(linspace(0,1.1,10))');
            set(hPressure(1),'LineWidth',3,'Color','b');
            set(hPressure(2),'LineWidth',2,'Color','k');   
        elseif length(Parameter) == 1 && Case == 2
            % Expansion fan
            % Expansion fans are gradual, so the solution should be 
            % a smooth curve
            xnew = [x;(x(end)+0.001:0.001:3)'];
            Z = length(xnew);
            LSpline = 400;
            S1 = [x(end-10:end);xnew(end-LSpline:end)];
            S2 = [PresRatio2(end-10:end);NPR*ones(LSpline+1,1)];
            A = nan(Z-L,1);
            for z = 1:length(A)
                A(z) = spline(S1,S2,xnew(L+z));
            end
            PresRatio2new = [PresRatio2;A];
            hPressure = plot(hpPressure,xnew,PresRatio2new,x(end)*ones(10,1),(linspace(0,1.1,10))');
            set(hPressure(1),'LineWidth',3,'Color','b');
            set(hPressure(2),'LineWidth',2,'Color','k');
        elseif length(Parameter) ~= 1 && Case == 0
            % Subsonic but not choked
            % Feed in the Mach number distribution 
            PresRatioSS = nan(L,1);
            for s = 1:L
                PresRatioSS(s) = PR(Parameter(s));
            end
            xnew = [x;(x(end)+0.001:0.001:3)'];
            Z = length(xnew);
            PresRatioSS = [PresRatioSS;PresRatioSS(end)*ones(Z-L,1)];
            hPressure = plot(hpPressure,xnew,PresRatioSS,x(end)*ones(10,1),(linspace(0,1.1,10))');
            set(hPressure(1),'LineWidth',3,'Color','b');
            set(hPressure(2),'LineWidth',2,'Color','k');
        elseif Case == -1
            % No flow
            xnew = [x;(x(end)+0.001:0.001:3)'];
            Z = length(xnew);
            PressureZero = ones(Z,1);
            hPressure = plot(hpPressure,xnew,PressureZero,x(end)*ones(10,1),(linspace(0,1.1,10))');
            set(hPressure(1),'LineWidth',3,'Color','b');
            set(hPressure(2),'LineWidth',2,'Color','k');
        elseif Parameter == true && Case == 0
            % Plot supersonic solution
            xnew = [x;(x(end)+0.001:0.001:3)'];
            Z = length(xnew);
            PresRatio2new = [PresRatio2;PresRatio2(end)*ones(Z-L,1)];
            hPressure = plot(hpPressure,xnew,PresRatio2new,x(end)*ones(10,1),(linspace(0,1.1,10))');
            set(hPressure(1),'LineWidth',3,'Color','b');
            set(hPressure(2),'LineWidth',2,'Color','k');     
        elseif Parameter == false && Case == 0
            % Plot subsonic solution 
            xnew = [x;(x(end)+0.001:0.001:3)'];
            Z = length(xnew);
            PresRationew = [PresRatio;PresRatio(end)*ones(Z-L,1)];
            hPressure = plot(hpPressure,xnew,PresRationew,x(end)*ones(10,1),(linspace(0,1.1,10))');
            set(hPressure(1),'LineWidth',3,'Color','b');
            set(hPressure(2),'LineWidth',2,'Color','k');               
        end
        title(hpPressure,'Pressure ratio plot');
        legend(hpPressure,'Current Solution','Nozzle Exit','Location','SouthWest');
        set(get(hpPressure,'XLabel'),'String','Horizontal Coordinate');
        set(get(hpPressure,'YLabel'),'String','Pressure Ratio');
        set(hpPressure,'XLim',[x(1) 1.5*x(end)],'YLim',[0 1.1]);
        grid(hpPressure);
    end

% Connected with GUI
p0 = 1500e6; % Some value % Only ratios are computed, so the numerical 
% value is not significant
tol = 1e-3;
% Analysis
if NPR > Subsonic && NPR ~= 1 && AR >= 1
    % Subsonic flow everywhere. 
    % Find the new sonic throat area based on the given value of NPR
    % DONE
    message = uicontrol('Style','text','String','Subsonic flow everywhere.','Position',[20,330,250,20]);
    set(message,'FontSize',11);
    NewSRA = FindSRA(NPR);
    MachSS = nan(L,1);
    for t = 1:L
        f = @(M) (1/M^2)*((2/(g+1))*(1+((g-1)/2)*M^2))^((g+1)/(g-1))... 
        - (radius(t)/NewSRA)^2;
        MachSS(t) = bisect(f,0.001,1,100);
    end
    PlotNozzle(0);
    PlotMach(MachSS,NPR,0);
    PlotPressure(MachSS,NPR,0);
elseif NPR == 1 && AR >= 1
    % No flow
    PlotNozzle(0)
    message = uicontrol('Style','text','String','No flow.','Position',[20,330,250,20]);
    set(message,'FontSize',11);
    PlotMach(2,NPR,-1);
    PlotPressure(2,NPR,-1);
    % Mach and pressure plots should be constant
elseif abs(NPR-Subsonic) < tol && AR >= 1
    % Subsonic solution 
    message = uicontrol('Style','text','String','Subsonic and choked.','Position',[20,330,250,20]);
    set(message,'FontSize',11);
    PlotNozzle(0);
    PlotMach(0,NPR,0);
    PlotPressure(0,NPR,0);
elseif NPR < Subsonic && NPR > Supersonic && AR > 1
    % Normal shock or oblique shock
    % Finding shock location 
    Mexit = sqrt((-1/(g-1))+sqrt((1/((g-1)^2))+(2/(g-1))*(2/(g+1))^((g+1)/(g-1))*(1/NPR)^2*(1/AreaRatio)^2));
    Pexit = NPR*p0; 
    Pstagexit = Pexit*(1+((g-1)/2)*Mexit^2)^(g/(g-1));
    Ratio = (Pstagexit/p0);
    PressureRatio = @(M) 1+(2*g/(g+1))*(M^2-1);
    MachRatioSQ = @(M) (M^2+(2/(g-1)))/((2*g/(g-1))*M^2-1);
    S = @(M) (PressureRatio(M))*(((1+((g-1)/2)*MachRatioSQ(M)))/(1+((g-1)/2)*M^2))^(g/(g-1)) - Ratio;
    MachCritical = bisect(S,0.9,5,100);
    tol = 1e-2;
    Index = nan;
    for q = 1:L
        if abs(Mach(q)-MachCritical) < tol
            Index = q;
            break
        end
    end
    if isnan(Index) == 1 && MachCritical > Mach(end) || Mexit > 1
    % Oblique shock
    message = uicontrol('Style','text','String','Oblique shock outside nozzle.','Position',[20,330,250,20]);
    set(message,'FontSize',11);
    % Oblique shock analysis
    % Need to find shock angle, final pressure ratio, and final Mach number
    FindPR = @(B) 1+(2*g/(g+1))*(Mach(end)^2*sin(B)^2-1) - NPR*(1/PresRatio2(end));
    beta = bisect(FindPR,0,pi/2,100);
    % Found shock angle in radians.
    % beta -- shock angle
    % Get theta, the angle the flow turns through. 
    theta = 2*cot(beta)*((Mach(end)^2*sin(beta)^2-1)/(Mach(end)^2*(g+cos(2*beta))+2));
    theta = tan(theta);
    % Found theta
    % Finally, get the exit Mach number
    MexitOS = (1/sin(beta-theta))*sqrt((1+((g-1)/2)*Mach(end)^2*sin(beta)^2)/(g*Mach(end)^2*sin(beta)^2-((g-1)/2)));
    % For PlotNozzle, feed in the shock angle, beta.
    % For PlotMach, feed in MexitOS
    % For PlotPressure, feed in NPR
    PlotNozzle(3,beta);
    PlotPressure(PresRatio2,NPR,2);
    PlotMach(Mach,MexitOS,2);
    else
        % Normal shock location found
    message = uicontrol('Style','text','String','Normal shock inside nozzle.','Position',[20,330,250,20]);
    set(message,'FontSize',11);
       ShockLocation = x(Index);
       PlotNozzle(1,ShockLocation,Index);
       MachNS = nan(L,1);
          MachNS(1:Index) = Mach(1:Index);
        % Based on Mexit, find the new sonic reference area
        ARnew = sqrt((1/(Mexit^2))*((2/(g+1))*(1+((g-1)/2)*Mexit^2))^((g+1)/(g-1)));
          SRA = rmin*AreaRatio/ARnew;
          for u = Index+1:L
              f = @(M) (1/M^2)*((2/(g+1))*(1+((g-1)/2)*M^2))^((g+1)/(g-1))... 
               - (radius(u)/SRA)^2;
              MachNS(u) = bisect(f,0.001,1,100);
          end
            PresRatioNS = nan(L,1);
            for w = 1:L
                if w > Index
                    PresRatioNS(w) = PR(MachNS(w))*Ratio;
                else
                    PresRatioNS(w) = PR(MachNS(w));
                end
            end
        PlotMach(MachNS,NPR,1);
        PlotPressure(PresRatioNS,NPR,1);
    end

elseif abs(NPR-Supersonic) < tol && AR > 1
    % Perfectly expanded
    message = uicontrol('Style','text','String','Perfectly expanded jet.','Position',[20,330,250,20]);
    set(message,'FontSize',11);
    PlotNozzle(0);
    PlotPressure(1,NPR,0);
    PlotMach(1,NPR,0);
else
    % Underexpanded jet
    % Prandtl-Meyer expansion fan
    message = uicontrol('Style','text','String','Underexpanded jet.','Position',[20,330,250,20]);
    set(message,'FontSize',11);
    Mout = PlotMach(true,NPR,2);
    PMAngle = (Prandtl_Meyer(Mout)-Prandtl_Meyer(Mach(end)));
    PlotPressure(true,NPR,2);
    PlotNozzle(2,PMAngle,Mach(end),Mout);
end
end
end