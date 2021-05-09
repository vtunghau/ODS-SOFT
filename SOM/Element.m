%------------------------------------------------------------------------------
%
% Element: Computes the osculating Keplerian elements from the satellite state
%          vector for elliptic orbits
%
% Inputs:
%   gm        Gravitational coefficient
%             (gravitational constant * mass of central body)
%   y         State vector (x,y,z,vx,vy,vz)
% 
% Outputs:
%             Keplerian elements (a,e,i,Omega,omega,M) with
%               a      Semimajor axis 
%               e      Eccentricity 
%               i      Inclination [rad]
%               Omega  Longitude of the ascending node [rad]
%               omega  Argument of pericenter  [rad]
%               M      Mean anomaly  [rad]
%
% Notes:
%   The state vector and gm must be given in consistent units, 
%   e.g. [m], [m/s] and [m^3/s^2]. The resulting unit of the semimajor
%   axis is implied by the unity of y, e.g. [m].
%
%   The function cannot be used with state vectors describing a circular
%   or non-inclined orbit.
%
% Last modified:   2015/08/12   M. Mahooti
%
%------------------------------------------------------------------------------
  function Kep = Element(gm, y)
  
  pi2 = 2*pi;
  
  r = y(1:3);                                  		 % Position
  v = y(4:6);                                  		 % Velocity
  
  h = cross(r,v);                                    % Areal velocity
  H = norm(h);
  
  Omega = atan2 ( h(1), -h(2) );                     % Long. ascend. node 
  Omega = mod(Omega,pi2);
  i     = atan2 ( sqrt(h(1)*h(1)+h(2)*h(2)), h(3) ); % Inclination        
  u     = atan2 ( r(3)*H, -r(1)*h(2)+r(2)*h(1) );    % Arg. of latitude   
  
  R  = norm(r);                                      % Distance           
  
  a = 1/(2/R-dot(v,v)/gm);                           % Semi-major axis    
  
  eCosE = 1-R/a;                                     % e*cos(E)           
  eSinE = dot(r,v)/sqrt(gm*a);                       % e*sin(E)           
  
  e2 = eCosE*eCosE +eSinE*eSinE;
  e  = sqrt(e2);                                     % Eccentricity 
  E  = atan2(eSinE,eCosE);                           % Eccentric anomaly  
  
  M  = mod(E-eSinE,pi2);                             % Mean anomaly
  
  nu = atan2(sqrt(1-e2)*eSinE, eCosE-e2);            % True anomaly
  
  omega = mod(u-nu,pi2);                             % Arg. of perihelion 
  
  Kep = [ a, e, i, Omega, omega, M ];
  
  