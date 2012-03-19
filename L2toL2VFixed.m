function [isfeasp,sOut] = L2toL2VFixed(f, h, x, W, V, R2, Gamma, sBasis, SOSoptions, SupplyRate)
%   function [isfeasp,sOut] = L2toL2VFixed(f, h, x, W, V, R2, Gamma, sBasis, SOSoptions, SupplyRate)
%
%   DESCRIPTION
%       This function solves the coupled SOS feasibility problems:
%          1. -Vdot*f +W^T*W - (1/Gamma)*h^T*h is SOS for {x: V(x) \leq R2)
%       Everything else is fixed except the SOS multipler for V(x) \leq R2
%
%   INPUTS expected data types indicated by [DATA TYPES]
%       f: state dynamics [polynomial]
%       h: output [polynomial]
%       x: state [polynomial]
%       W: input [polynomial]
%       V: storage function for differential inequality [polynomial]
%       R2: bound on level set of V in {x: V(x) \leq R2} [double]
%       Gamma: bound on gain [double]
%       sBasis: polynomial basis for decision variable s [polynomial]
%       SOSoptions: options, (optional) [sosoptions]
%       SupplyRate: (optional) If the system is dissipative wrt supply rates, they are
%          listed as a cell of polynomials such that:
%          \int_0^T SupplyRate{i} dt \geq 0 for all T>0 , for i = 1 to N%           
%           [1xN cell of polynomials]
%
%   OUTPUTS:
%       isfeasp: feasibility results [boolean]
%       sOut: sos-multiplier for set containment {x: V(x) \leq R2} [polynomial]
if(nargin == nargin(@L2toL2VFixed)-1)
    SupplyRate = [];
end

if(nargin == nargin(@L2toL2VFixed)-2)
    SOSoptions = sosoptions();
end

% initilaization
sOut = polynomial(0);
s = polydecvar('c',sBasis,'mat');
Vdot = jacobian(V,x)*f;
RobustTerm = polynomial(0);

if(~isempty(SupplyRate))
    %if a supply rates are included, add each rate to the differential
    %inequality, as well as sos multipliers
    for i=1:length(SupplyRate)
        sR{i} = polydecvar(sprintf('r%d', i), monomials([x; W],0:0), 'mat');
        RobustTerm = RobustTerm + sR{i}*SupplyRate{i};
        sosineq{i} = sR{i};
    end    
end
%differential inequality
sosineq{end+1} = -((R2-V)*s + Vdot - W'*W + h'*h/Gamma^2) -RobustTerm;
% V(x) <R^2 set containment sos multiplier
sosineq{end+1} = s;

%SOS optimization
[info,dopt] = sosopt(sosineq,[x;W], [], SOSoptions);
isfeasp = info.feas;
if isfeasp == 1
    sOut = subs(s,dopt);
end