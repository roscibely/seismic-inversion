function [I_wagn, s_wagn]=invwagn(s,t,fb,fa, w, a_o, ini_i, iter_max)

%[irc,dp] = wagn(s,w,a_o,iter_max); Quasi Newton
[irc,dp] = gradd(s,w,a_o,iter_max,10e-6);
%a_o = 0.08;
%[irc,dp] = wagn(s,w,a_o,iter_max);
delf=min([5 fb/5]);
I1= exp(2*cumsum(irc));
I1=I1-mean(I1);
p=polyfit(t,ini_i,1);
implog=ini_i-polyval(p,t);
implog=pad(implog,I1);%zero pad
t2=xcoord(t(1),t(2)-t(1),implog);
%
imp=mergetrcs(implog,I1,t2,fb,delf,fa);
imp=imp(1:length(s));
%
imp=imp+polyval(p,t);
I_wagn = imp;
s_wagn = dp;
end 