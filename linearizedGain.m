function [gL,Vout] = linearizedGain(f,h,x,w,l1, SOSOptions)
% Copywright (c) Erin Summers 2012
%
%function [gL,Vout] = linearizedGain(f,h,x,w,l1, SOSOptions)
%   DESCRIPTION
%       This function finds an L2 gain for which the linearized dynamics of
%       a system, which is valid for a small ball around the origin
%
%   INPUTS expected data types indicated by [DATA TYPES]
%       f: state dynamics [polynomial]
%       h: output [polynomial]
%       x: state [polynomial]
%       w: input [polynomial]
%       l1: some small polynomial, usually of the form eps*x'x
%           [polynomial], ensures V is decresent and positive definite
%       SOSoptions: options, (optional) [sosoptions]
%
%   OUTPUTS:
%       gL: L2 gain returned [double]
%       Vout: storage function that certifies L2 gain

if(nargin<6)
    SOSOptions = sosoptions();
end

%initialization
gL = inf;
Vout = polynomial(0);
V = polydecvar('c',monomials(x,2),'vec');

%get linearization of system
fLin = cleanpoly(f,[],1);
hLin = cleanpoly(h,[],1);
A = jacobian(fLin,x);
B = jacobian(fLin,w);
C = jacobian(hLin,x);

%bisection bounds
gsUp = 1e3;
gsLow = 0;
gs = gsUp;

%positive semidefinite constraint
sosineq{1} = V-l1;

%differential inequality constraint
sosineq{2} = gs*w'*w-x'*C'*C*x - jacobian(V,x)*(A*x+B*w);

%SOS optimization
[info,dopt]=sosopt(sosineq,[x;w], [],SOSOptions);
isfeasp = info.feas;
%if first gain gL is feasbile, return, else bisect on gL
if isfeasp
    gL = sqrt(gs);
    Vout = subs(V,dopt);
end

%bisect on gL
while ~isfeasp
    
    gsLow = gsUp;
    gsUp = 2*gsUp;
    gs = gsUp;
    
    
    V = polydecvar('c',monomials(x,2),'vec');
    sosineq{1} = V-l1;
    sosineq{2} = gs*w'*w-x'*C'*C*x - jacobian(V,x)*(A*x+B*w);
    [info,dopt] = sosopt(sosineq,[x;w],[],SOSOptions);
    isfeasp = info.feas;
    if isfeasp
        gL = sqrt(gs);
        Vout = subs(V,dopt);
    end
end

while gsUp - gsLow > 1e-3
    gs = 0.5*(gsUp+gsLow);
        V = polydecvar('c',monomials(x,2),'vec');
    sosineq{1} = V-l1;
    sosineq{2} = gs*w'*w-x'*C'*C*x - jacobian(V,x)*(A*x+B*w);
    [info,dopt] = sosopt(sosineq,[x;w]);
    isfeasp = info.feas;
    if isfeasp
        gsUp = gs;
        gL = sqrt(gs);
        Vout = subs(V,dopt);
    else
        gsLow = gs;
    end
end

