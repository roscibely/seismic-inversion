function [I_recursiva, r_recursiva]=invrecursiva(s,ini_i,t,fb,fa, w)

delf=min([5 fb/5]);
%irc=blint(s,t,fb,fa); %filtro gaussiano
[f,o] = spiking(s,20,0.9); %filtro Wiener
irc = conv2(f,s);

I1= exp(2*irc );
I1=I1-mean(I1);
p=polyfit(t,ini_i,1);
implog=ini_i-polyval(p,t);
implog=pad(implog,I1);%zero pad
t2=xcoord(t(1),t(2)-t(1),implog);
imp=mergetrcs(implog,I1,t2,fb,delf,fa);
imp=imp(1:length(s));
imp=imp+polyval(p,t);
I_recursiva = imp;
r_recursiva = irc;
%s_recursiva = conv2(r_recursiva, w);
end