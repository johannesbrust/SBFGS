clc;
clear;
clear all;
clear all hidden;

num_result = 3;

path_full = '../RESULT_Paper/OK/SBFGS_FullRank_addNLP_2var/';
path_low = '../RESULT_Paper/OK/SBFGS_LowRank_addNLP_2var/'; 

path = path_full;

file1 = strcat(path,'stats_hessian_descent.mat');
file2 = strcat(path,'stats_bfgs.mat');
file3 = strcat(path,'stats_sbfgs_kplus1_inertia.mat');
file4 = strcat(path,'stats_sbfgs_kplus1_descent.mat');

file = { file1, file2, file3};

maxLen = 0;
for i=1:num_result
    load(char(file(i)));
    solAll{i} = sol;
    if maxLen<length(sol)
      maxLen = length(sol);
      baseSol = sol;
    end
end

idx = 0;
fprintf('----> Running %d cases:\n',length(sol))   ;
fprintf('                          #Hes      |        #BFGS     |        #SBFGS_P  | \n')  ; 
fprintf('      Prob         #St    #It   #Fun| #St    #It   #Fun| #St    #It   #Fun| \n')  ; 


sol_ok = zeros(num_result,1);
sol_max = zeros(num_result,1);
sol_id = ones(num_result,1);

sol_stat  = zeros(maxLen,num_result); 
sol_iter  = zeros(maxLen,num_result); 
sol_efunc = zeros(maxLen,num_result); 
    
for idx = 1:maxLen
    for nSet = 1:num_result
       prob_name = {solAll{nSet}.name};
       prob_iter = [solAll{nSet}.nIter];
       prob_efun = [solAll{nSet}.zoomLS];
       prob_flag = {solAll{nSet}.flag};   
       if strcmp(baseSol(idx).name,prob_name(sol_id(nSet)))
            sol_iter(idx,nSet)  = prob_iter(sol_id(nSet));
            sol_efunc(idx,nSet) = prob_efun(sol_id(nSet));
            if strcmp(prob_flag(sol_id(nSet)),'OK  ')
                sol_ok(nSet) = sol_ok(nSet)+1;
                sol_stat(idx,nSet) = 1;
            elseif strcmp(prob_flag(sol_id(nSet)),'MAX ')    
                sol_max(nSet) = sol_max(nSet)+1;
                sol_stat(idx,nSet) = 2;
            end
            sol_id(nSet) = sol_id(nSet)+1;
        end    
    end
    
    fprintf('%16s%3d%6d%6d|%3d%6d%6d|%3d%6d%6d \n',baseSol(idx).name,...
    sol_stat(idx,1), sol_iter(idx,1),sol_efunc(idx,1),sol_stat(idx,2), sol_iter(idx,2),sol_efunc(idx,2),...
    sol_stat(idx,3), sol_iter(idx,3),sol_efunc(idx,3));
end

fprintf('\n');

fprintf('                OK      MAX      FAIL   Iter       Func     F/I\n');
fprintf('Total:           %d\n', maxLen);
fprintf('       Hes       %d      %d       %d     %d       %d     %2.2f\n', sol_ok(1), sol_max(1), maxLen-sol_ok(1)-sol_max(1), sum(sol_iter(:,1)), sum(sol_efunc(:,1)),sum(sol_efunc(:,1))/sum(sol_iter(:,1)));
fprintf('       Bfgs      %d      %d       %d     %d       %d     %2.2f\n', sol_ok(2), sol_max(2), maxLen-sol_ok(2)-sol_max(2), sum(sol_iter(:,2)), sum(sol_efunc(:,2)),sum(sol_efunc(:,2))/sum(sol_iter(:,2)));
fprintf('       Sbfgs_P   %d      %d       %d     %d       %d     %2.2f\n', sol_ok(3), sol_max(3), maxLen-sol_ok(3)-sol_max(3), sum(sol_iter(:,3)), sum(sol_efunc(:,3)),sum(sol_efunc(:,3))/sum(sol_iter(:,3)));


for ii=1:num_result; ki(ii)=sum(sol_iter(sol_stat(:,ii)==1,ii)); end
for ii=1:num_result; kf(ii)=sum(sol_efunc(sol_stat(:,ii)==1,ii)); end

lessi = zeros(num_result,1);
lessf = zeros(num_result,1);
totalp = zeros(num_result,1);
for nSet = 3:num_result 
    for idx = 1:maxLen
		if sol_stat(idx,nSet)==1 && ( sol_stat(idx,2)==1)
			totalp(nSet) = totalp(nSet)+1;
	   		if ( sol_iter(idx,nSet) <= sol_iter(idx,2) )
	   	  		lessi(nSet) = lessi(nSet)+1;
			end
	   		if ( sol_efunc(idx,nSet) <= sol_efunc(idx,2) )
	   	  		lessf(nSet) = lessf(nSet)+1;
			end
	   	end
	end
end


fprintf('\n only for OK prob: \n');
fprintf('                  OK      Iter       Func     F/I      LessIter   LessEval\n');
fprintf('       Hes        %d      %d       %d     %2.2f\n', sol_ok(1), ki(1), kf(1),kf(1)/ki(1));
fprintf('       Bfgs       %d      %d       %d     %2.2f\n', sol_ok(2), ki(2), kf(2),kf(2)/ki(2));
fprintf('       Sbfgs_P    %d      %d       %d     %2.2f      %d/%d      %d/%d\n', sol_ok(3), ki(3), kf(3),kf(3)/ki(3),lessi(3), totalp(3),lessf(3), totalp(3));



%%% generate figure
num_prob = maxLen;

IterT = zeros(num_prob,num_result);
EfunT = zeros(num_prob,num_result);

for i=1:num_prob
    for j=1:num_result
        if sol_stat(i,j) ~= 1
           IterT(i,j) = nan;
           EfunT(i,j) = nan;           
        else
           IterT(i,j) = sol_iter(i,j);
           EfunT(i,j) = sol_efunc(i,j);           
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

xlabel(' log_2(\tau) - Iterations','FontSize',14)
ylabel('% Problems','FontSize',14)
axis([0 7.5 0 100])
grid on
legend('Hes','Bfgs', 'Sbfgs\_P',4)
legend('boxoff')
set(gca,'FontSize',14)

if strcmp(path, path_full)
    print -depsc iter_fullrank_nlp2var_kplus1.eps
else
    print -depsc iter_lowrank_nlp2var_kplus1.eps
end



[np,ns,r]=perf(EfunT,1);
figure(2)
[xs,ys] = stairs(r(:,1),[1:np]/np);
plot(xs,100*ys,'b','MarkerSize',3);
hold on;
[xs,ys] = stairs(r(:,2),[1:np]/np);
plot(xs,100*ys,'r','MarkerSize',3);
[xs,ys] = stairs(r(:,3),[1:np]/np);
plot(xs,100*ys,'g','MarkerSize',3);

xlabel(' log_2(\tau) -  Function Evaluations','FontSize',14)
ylabel('% Problems','FontSize',14)
axis([0 7.5 0 100])
grid on
legend('Hes','Bfgs', 'Sbfgs\_P',4)
legend('boxoff')
set(gca,'FontSize',14)

if strcmp(path, path_full)
    print -depsc efunc_fullrank_nlp2var_kplus1.eps
else
    print -depsc efunc_lowrank_nlp2var_kplus1.eps
end


        