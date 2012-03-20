function [Vout] = linearizedGainGivenGamma(f,h,x,w,l1, gL,SOSoptions)
%function [Vout] = linearizedGainGivenGamma(f,h,x,w,l1, gL,SOSoptions)
%   DESCRIPTION
%       This function solves the coupled feasibility problems:
%           1. -Vdot*fLin + gL^2w'*w - hLin'*hLin
%           2. V is positive definite
%           where V is free, and everything else is fixed
%
%   INPUTS expected data types indicated by [DATA TYPES]
%       f: state dynamics [polynomial]
%       h: output [polynomial]
%       x: state [polynomial]
%       w: input [polynomial]
%       l1: some small polynomial, usually of the form eps*x'x
%           [polynomial], ensures V is decresent and positive definite
%       gL: proposed L2 gain [double]
%       SOSoptions: options, (optional) [sosoptions]
%
%   OUTPUTS:
%       Vout: storage function that certifies L2 gain

if(nargin<7)
    SOSoptions = sosoptions();
end

%initialization, gain is squared in SOS problem
Vout = polynomial(0);
V = polydecvar('c',monomials(x,2),'vec');
gs = gL^2;

%find linearization
fLin = cleanpoly(f,[],1);
hLin = cleanpoly(h,[],1);
A = jacobian(fLin,x);
B = jacobian(fLin,w);
C = jacobian(hLin,x);


%positive definite constraint
sosineq{1} = V-l1;

%differential inequaliy constraint
sosineq{2} = gs*w'*w-x'*C'*C*x - jacobian(V,x)*(A*x+B*w);

%SOS optimization
[info,dopt]=sosopt(sosineq,[x;w],[],SOSoptions);
isfeasp = info.feas;

%if result is feasible, return the result
if isfeasp
    Vout = subs(V,dopt);
end
