function argOut = myBisection( xL_Lim, xL, xH, xH_Lim, fh, fargin, bisArg, BisTol, allowMAX)
%function argOut = myBisection( xL_Lim, xL, xH, xH_Lim, fh, fargin, bisArg, BisTol, allowMAX)
% DESCRIPTION
%   This function is performs a bisection by choosing useful bisection
%   bounds and limits of the bounds. In order for the bisection to run, the
%   lower bound (xL) must be feasible for the function fh, and the upper
%   bound (xH) must be infeasible.  If the inital arguments passed into
%   this function don't meet this criteria, then xL and xH are pushed
%   towards xL_Lim and xH_Lim, respectively, until the pass/fail criteria
%   is met (or else the function terminates).
% INPUTS
%   xL_Lim : absolute limit of lower bound. If the bisection bound is
%       pushed below this point, the function terminates.
%   xL: lower bound of bisection. The sos program in fh must be feasible
%       when evaluated at xL.
%   xH: upper bound of bisection. The sos program in fh must be infeasible
%       when evaluated at xH.
%   xH_Lim: absolute limit of upper bound. If the bisection bound is
%       pushed above this point, the function terminates.
%   fh: function handle of an sos optimiaztion program.This function
%       determines if the optimization evaluated with the bisection
%       varaible is feasible or infeasible.  It is assumed that the 1st
%       output argument of fh is a feasibility flag.  The remaining output
%       arguments are polynomial objects returned from the sos
%       optimization.  myBisection is flexible and adapts to a varying
%       number of output arguments.
%   fargin: input arguments to function handle fh
%   bisArg: the position of the input argument to fh for which the
%       bisection will be done.  For example, if you which to bisect on the
%       5th input argument of fh, then bisArg = 5;
%   BisTol: relative bisection tolerance
%   allowMAX: optional, allow analysis to continue even if upper bound is
%       viloated
% Outputs:
%   argOut: the first argument of argOut is the final value of the
%   bisection variable.  The remaining arguments are the output(s) from the
%   sos program in fh

if(nargin < 9)
    allowMAX = false;
end

if(xL_Lim ==0)
    xL_Lim = eps;    
end
if(xL ==0)
    xL = eps;
end

for i = 1:nargout(fh)
    evalc(['argOut{i} = 0']);
end
% Sanity Check
if ~(xL_Lim <=xL && xL <=xH && xH<=xH_Lim)
    fprintf('\n One of the inequalities was viloated: xL_Lim <=xL <= xH <= xH_Lim. Fix your inputs!\n\n')
    argOut{1} = 0;
    for i = 2:nargout(fh)
         argOut{i} = 0;
    end

    return;
end

%Establish xL
xLtry = xL;
lowFactor = 2;
highFactor = 2;
xH_estab = 0;

fprintf('Establishing Lower Bisection Bound\n')
%fprintf('xL_Lim \t xLtry \t xH \t xH_Lim\n \t feasible?')
xL_Limfeas = 0;
while(xLtry >= xL_Lim)
    fargin{bisArg} = xLtry;
    isfeasp = fh(fargin{:});
    if(isfeasp)
        xL = xLtry;
        xL_Limfeas = 1;
        %fprintf('%6.4f \t %6.4f \t %6.4f \t %6.4f \n', xL_Lim, xLtry, xH, xH_Lim);
        break;
    else
        xH_Lim = xLtry;
        xH = xLtry;
        xHtry = xLtry;
        xLtry = xLtry/lowFactor;
        xH_estab = 1;
    end
    %fprintf('%6.4f \t %6.4f \t %6.4f \t %6.4f \n', xL_Lim, xLtry, xH, xH_Lim);
end

%try xL_Lim if you haven't found anything feasible, since it should be feasible
if(~xL_Limfeas)
    xtry = xL_Lim;
    fargin{bisArg} = xtry;
    isfeasp = fh(fargin{:});
    if(isfeasp)
        xL = xtry;
    else 
        fprintf('Limit of Lower Bisection Bound is not feasible! Lower it!\n')        
        return;
    end
end

% if(xLtry<xL_Lim)
%     fprintf('Limit of Lower Bisection Bound violated\n')
%     return;
% end

%Establish xH
fprintf('Establishing Upper Bisection Bound\n')
%fprintf('xH \t xH_Lim \t xHtry \n')
if(~xH_estab)
    xHtry = xH;
    while(xHtry <=xH_Lim)
        fargin{bisArg} = xHtry;
        isfeasp = fh(fargin{:});
        if(isfeasp)
            %fprintf('%6.4f \t %6.4f \t %6.4f\n', xH, xH_Lim, xHtry)
            xHtry = xHtry*highFactor;
            %xL = xHtry; %es 2-4-11
        else
            %fprintf('%6.4f \t %6.4f \t %6.4f\n', xH, xH_Lim, xHtry)
            xH = xHtry;
            break;
        end
    end
end

if(xHtry>xH_Lim)
    if(~allowMAX)
        fprintf('Limit of Upper Bisection Bound violated\n')
    else
        xtry = xH_Lim;
        fargin{bisArg} = xtry;
        evalFun = '[isfeasp';
        for i = 2:nargout(fh)
            evalFun  = [evalFun, ' arg', num2str(i)];
        end
        evalFun = [evalFun '] = fh(fargin{:})'];
        evalc(evalFun);
        fprintf('%6.4f \t %6.4f \t  %6.4f \t %i\n', xtry, xL, xH, isfeasp)
        argOut{1} = xtry;
        for i = 2:nargout(fh)
            evalc(['argOut{i} = arg' num2str(i)]);
        end
    end
    return;
end

% Run Bisection
fprintf('Running Bisection\n')
fprintf('xtry \t xLow \t \t xHigh \t \t isfeas\n')
xtry = (xL+xH)/2;
CurTol = 1;
argOut{1} = 0;
never_feas = 1;
for i = 2:nargout(fh)
    evalc(['argOut{i} = 0']);
end


CurTol = 1;
firstiter = 1;

while(CurTol> BisTol)
    fargin{bisArg} = xtry;    
    evalFun = '[isfeasp';
    for i = 2:nargout(fh)
        evalFun  = [evalFun, ' arg', num2str(i)];
    end
    evalFun = [evalFun '] = fh(fargin{:})'];
    evalc(evalFun);
    
    fprintf('%6.4f \t %6.4f \t  %6.4f \t %i\n', xtry, xL, xH, isfeasp)
    if ~isfeasp
        xH = xtry;
    else
        never_feas = 0;
        xL = xtry;
        argOut{1} = xtry;
        for i = 2:nargout(fh)
            evalc(['argOut{i} = arg' num2str(i)]);
        end
    end
    
        xOld = xtry;
    xtry = 0.5*(xL + xH);
    CurTol = (xH-xtry)/xtry;

end
%in the case where none of the points tested are feasible, the xL
%should at least be feasible
if(never_feas)
    xtry = xL;
    fargin{bisArg} = xtry;
    
    evalFun = '[isfeasp';
    for i = 2:nargout(fh)
        evalFun  = [evalFun, ' arg', num2str(i)];
    end
    evalFun = [evalFun '] = fh(fargin{:})'];
    evalc(evalFun);
    fprintf('%6.4f \t %6.4f \t  %6.4f \t %i\n', xtry, xL, xH, isfeasp)
    argOut{1} = xtry;
    for i = 2:nargout(fh)
        evalc(['argOut{i} = arg' num2str(i)]);
    end
end