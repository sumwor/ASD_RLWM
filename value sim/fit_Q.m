function [qpar, negloglike]=fit_Q(stats,fit_fun,initpar)    
% % fit_Q %
%PURPOSE:   Fit the choice behavior to a q-learning model based on reward prediction error
%           using maximum likelihood estimate
%AUTHORS:   AC Kwan 170518
%
%INPUT ARGUMENTS
%   stats:      stats of the game
%   fit_fun:    the function to fit, e.g., Q_RPEfun
%   initpar:    initial values for the parameters
%   x:          which player? (1 or 2)
%
%OUTPUT ARGUMENTS
%   qpar:       extracted parameters (alpha and beta)
%   negloglike: negative log likelihood

%%
maxit=1e6;
maxeval=1e6;
op=optimset('fminsearch');
op.MaxIter=maxit;
op.MaxFunEvals=maxeval;

c = stats.c(~isnan(stats.c));  %take player x's choice history
r = stats.r(~isnan(stats.c));  %take player x's reward history
s = stats.s(~isnan(stats.c));

func_handle = str2func(fit_fun);
[qpar, negloglike, exitflag, output]=fminsearch(func_handle, initpar, [], [c r s]);
%[qpar, negloglike, exitflag, output]=fminsearch(@Q_RPEfun, initpar, [], [c r]);

if exitflag==0
    qpar=nan(size(qpar));   %did not converge to a solution, so return NaN
    negloglike=nan;
end

end
