function output = FindSRA(in)
global rmin AreaRatio g
     % Output is the new sonic reference area
     ExitMach = sqrt((2/(g-1))*((in)^(-(g-1)/g)-1));
     A_Astar = sqrt((1/(ExitMach^2))*((2/(g+1))*(1+((g-1)/2)*ExitMach^2))^((g+1)/(g-1)));
     % AR is A/A*
     output = rmin*AreaRatio/A_Astar;
end