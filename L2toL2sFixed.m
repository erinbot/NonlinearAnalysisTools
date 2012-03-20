function [isfeas,VOut] = L2toL2sFixed(f, h, x, W, s, R2, Gamma, L1, VBasis, SOSoptions, SupplyRate)
% Copywright (c) Erin Summers 2012
%
%   function [isfeas,VOut] = L2toL2sFixed(f, h, x, W, s, R2, Gamma, L1, VBasis, SOSoptions, SupplyRate)
%
%   DESCRIPTION
%       This function solves the coupled SOS feasibility problems:
%          1. -Vdot*f +W^T*W - (1/Gamma)*h^T*h is SOS for {x: V(x) \leq R2)
%          2. V is SOS
%       over choice of V. Everything else is fixed.
%
%   INPUTS expected data types indicated by [DATA TYPES]
%       f: state dynamics [polynomial]
%       h: output [polynomial]
%       x: state [polynomial]
%       W: input [polynomial]
%       s: sos-multiplier for set containment {x: V(x) \leq R2} [polynomial]
%       R2: bound on level set of V in {x: V(x) \leq R2} [double]
%       Gamma: bound on gain [double]
%       L1: some small polynomial, usually of the form eps*x'x
%           [polynomial], ensures V is decresent and positive definite
%       VBasis: polynomial basis for decision variable V [polynomial]
%       SOSoptions: options, (optional) [sosoptions]
%       SupplyRate: (optional) If the system is dissipative wrt supply rates, they are
%          listed as a cell of polynomials such that:
%          \int_0^T SupplyRate{i} dt \geq 0 for all T>0 , for i = 1 to N%           
%           [1xN cell of polynomials]
%
%   OUTPUTS:
%       isfeas: feasibility results [boolean]
%       VOut: result from SOS optimization, zero if infeasible [polynomial]

if(nargin == nargin(@L2toL2sFixed)-1)
    SupplyRate = [];
end

if(nargin == nargin(@L2toL2sFixed)-2)
    SOSoptions = sosoptions();
end

% initialization
VOut = polynomial(0);
V = polydecvar('c',VBasis,'mat');
Vdot = jacobian(V,x)*f;
sosineq = {};

% positive definite, decresent V constraint
sosineq{1} = V - L1;

RobustTerm = polynomial(0);
if(~isempty(SupplyRate))
    %if supply rate included, add sos multiplier and include it in the
    %differential inequality constraint;
     for i=1:length(SupplyRate)
        sR{i+1} = polydecvar(sprintf('r%d', i), monomials([x; W],0:0), 'mat');
        RobustTerm = RobustTerm + sR{i}*SupplyRate{i};
        sosineq{i+2} = sR{i};
    end
end
%differential inequality constraint
sosineq{end+1} = -((R2-V)*s + Vdot - W'*W + h'*h/Gamma^2) -RobustTerm;

% SOS optimization 
[info,dopt] = sosopt(sosineq,[x;W],[], SOSoptions);
isfeas = info.feas;

%if result is feasible, return it, else return zero polynomials
if isfeas == 1
    VOut = subs(V,dopt);
end