function [] = do_prof_NO_INF(maxNProb,nSet,iter,efunc,okflag,ctlflag)


num_result = nSet;
num_prob = maxNProb;

IterT = zeros(num_prob,num_result);
EfunT = zeros(num_prob,num_result);

for i=1:num_prob
    for j=1:num_result
        if okflag(i,j) ~= 1
           IterT(i,j) = nan;
           EfunT(i,j) = nan;           
        else
           IterT(i,j) = iter(i,j);
           EfunT(i,j) = efunc(i,j);           
        end            
    end
end


[np,ns,r]=perf(IterT,1);
figure(1)
[xs,ys] = stairs(r(:,1),[1:np]/np);
plot(xs,100*ys,'b','MarkerSize',3);
hold on;
[xs,ys] = stairs(r(:,2),[1:np]/np);
plot(xs,100*ys,'r','MarkerSize',3);
[xs,ys] = stairs(r(:,3),[1:np]/np);
plot(xs,100*ys,'g','MarkerSize',3);
[xs,ys] = stairs(r(:,4),[1:np]/np);
plot(xs,100*ys,'m','MarkerSize',3);
[xs,ys] = stairs(r(:,5),[1:np]/np);
plot(xs,100*ys,'c','MarkerSize',3);
[xs,ys] = stairs(r(:,6),[1:np]/np);
plot(xs,100*ys,'k','MarkerSize',3);

if nSet==7
    [xs,ys] = stairs(r(:,7),[1:np]/np);
    plot(xs,100*ys,'--','MarkerSize',3);
end
    
xlabel(' log_2(\tau) - Iterations','FontSize',14)
ylabel('% Problems','FontSize',14)
axis([0 7.5 0 100])
grid on
if nSet==6
	if ctlflag== 1
		legend('hes_i','hes_d','bfgs','sbfgs_1','sbfgs_1Li','sbfgs_1Ld',4)
	else
    	legend('hes_i','hes_d','bfgs','sbfgs_1','sbfgs_2i','sbfgs_2d',4)
	end
end
if nSet==7
    legend('hes_i','hes_d','bfgs','sbfgs_2i','sbfgs_2d','sbfgs_2Li','sbfgs_2Ld',4)
end
legend('boxoff')
set(gca,'FontSize',14)

if nSet==6 && ctlflag==1
% print -depsc iter_lowrank_nlp2var_k.eps
elseif ctlflag==1
 print -depsc iter_lowrank_nlp2var_kplus1.eps
end
%  print -depsc iter_fullrank_nlp2var.eps
%  print -depsc iter_fullrank_setRatio010.eps
%  print -depsc iter_fullrank_setRatio025.eps
%  print -depsc iter_fullrank_setRatio050.eps
%  print -depsc iter_fullrank_setRatio075.eps
%  print -depsc iter_fullrank_setRatio090.eps




[np,ns,r]=perf(EfunT,1);
figure(2)
[xs,ys] = stairs(r(:,1),[1:np]/np);
plot(xs,100*ys,'b','MarkerSize',3);
hold on;
[xs,ys] = stairs(r(:,2),[1:np]/np);
plot(xs,100*ys,'r','MarkerSize',3);
[xs,ys] = stairs(r(:,3),[1:np]/np);
plot(xs,100*ys,'g','MarkerSize',3);
[xs,ys] = stairs(r(:,4),[1:np]/np);
plot(xs,100*ys,'m','MarkerSize',3);
[xs,ys] = stairs(r(:,5),[1:np]/np);
plot(xs,100*ys,'c','MarkerSize',3);
[xs,ys] = stairs(r(:,6),[1:np]/np);
plot(xs,100*ys,'k','MarkerSize',3);

if nSet==7
    [xs,ys] = stairs(r(:,7),[1:np]/np);
    plot(xs,100*ys,'--','MarkerSize',3);
end
    
xlabel(' log_2(\tau) - Function Evaluations','FontSize',14)
ylabel('% Problems','FontSize',14)
axis([0 7.5 0 100])
grid on
if nSet==6
	if ctlflag== 1
		legend('hes_i','hes_d','bfgs','sbfgs_1','sbfgs_1Li','sbfgs_1Ld',4)
	else
    	legend('hes_i','hes_d','bfgs','sbfgs_1','sbfgs_2i','sbfgs_2d',4)
	end
end
if nSet==7
    legend('hes_i','hes_d','bfgs','sbfgs_2i','sbfgs_2d','sbfgs_2Li','sbfgs_2Ld',4)
end
legend('boxoff')
set(gca,'FontSize',14)


if nSet==6 && ctlflag==1
% print -depsc efun_lowrank_nlp2var_k.eps
elseif ctlflag==1
 print -depsc efun_lowrank_nlp2var_kplus1.eps
end
%  print -depsc efun_fullrank_nlp2var.eps
%  print -depsc efun_fullrank_setRatio010.eps
%  print -depsc efun_fullrank_setRatio025.eps
%  print -depsc efun_fullrank_setRatio050.eps
%  print -depsc efun_fullrank_setRatio075.eps
%  print -depsc efun_fullrank_setRatio090.eps

