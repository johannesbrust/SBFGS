clc;
clear;
clear all;
clear all hidden;

num_result = 10;

path = '../RESULT_Paper/OK/SBFGS_FullRank_addNLP_2var/';
path_low = '../RESULT_Paper/OK/SBFGS_LowRank_addNLP_2var/'; 
% path = '../RESULT_Paper/OK/setRatio_010/';
% path = '../RESULT_Paper/OK/setRatio_025/';
% path = '../RESULT_Paper/OK/setRatio_050/';
% path = '../RESULT_Paper/OK/setRatio_075/';
% path = '../RESULT_Paper/OK/setRatio_090/';

file1 = strcat(path,'stats_hessian_inertia.mat');
file2 = strcat(path,'stats_hessian_descent.mat');
file3 = strcat(path,'stats_bfgs.mat');
file4 = strcat(path,'stats_sbfgs_k.mat');
file5 = strcat(path,'stats_sbfgs_kplus1_inertia.mat');
file6 = strcat(path,'stats_sbfgs_kplus1_descent.mat');

file = { file1, file2, file3, file4, file5, file6};

if num_result == 10
	file7 = strcat(path_low,'stats_sbfgs_low_k_inertia.mat');
	file8 = strcat(path_low,'stats_sbfgs_low_k_descent.mat');
	file9 = strcat(path_low,'stats_sbfgs_low_kplus1_inertia.mat');
	file10 = strcat(path_low,'stats_sbfgs_low_kplus1_descent.mat');
	file = { file1, file2, file3, file4, file5, file6,file7, file8,file9, file10};
end



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
if num_result == 6
    fprintf('                          #Hes_I    |        #Hes_D      |         #BFGS      |        #SBFGS_1    |       #SBFGS_2i     |      #SBFGS_2d                \n')  ; 
    fprintf('      Prob         #St     #It  #Fun|   #St     #It  #Fun|   #St     #It  #Fun|   #St     #It  #Fun|   #St     #It  #Fun|   #St     #It  #Fun \n')  ; 
elseif  num_result == 10
    fprintf('                          #Hes_I    |        #Hes_D      |         #BFGS      |        #SBFGS_1    |      #SBFGS_1Li   |     #SBFGS_K_1Ld     |       #SBFGS_2i     |      #SBFGS_2d      |       #SBFGS_2Li     |      #SBFGS_2Ld               \n')  ; 
    fprintf('      Prob         #St     #It  #Fun|   #St     #It  #Fun|   #St     #It  #Fun|   #St     #It  #Fun|   #St     #It  #Fun|   #St     #It  #Fun|   #St     #It  #Fun|   #St     #It  #Fun|   #St     #It  #Fun|   #St     #It  #Fun \n')  ; 
end

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
 
if num_result == 6   
        fprintf('%16s%3d%6d%6d|%3d%6d%6d|%3d%6d%6d|%3d%6d%6d|%3d%6d%6d|%3d%6d%6d \n',baseSol(idx).name,...
        sol_stat(idx,1), sol_iter(idx,1),sol_efunc(idx,1),sol_stat(idx,2), sol_iter(idx,2),sol_efunc(idx,2),...
        sol_stat(idx,3), sol_iter(idx,3),sol_efunc(idx,3),sol_stat(idx,4), sol_iter(idx,4),sol_efunc(idx,4),...
        sol_stat(idx,5), sol_iter(idx,5),sol_efunc(idx,5),sol_stat(idx,6), sol_iter(idx,6),sol_efunc(idx,6));
elseif num_result == 10
        fprintf('%16s%3d%6d%6d|%3d%6d%6d|%3d%6d%6d|%3d%6d%6d|%3d%6d%6d|%3d%6d%6d|%3d%6d%6d|%3d%6d%6d|%3d%6d%6d|%3d%6d%6d \n',baseSol(idx).name,...
        sol_stat(idx,1), sol_iter(idx,1),sol_efunc(idx,1),sol_stat(idx,2), sol_iter(idx,2),sol_efunc(idx,2),...
        sol_stat(idx,3), sol_iter(idx,3),sol_efunc(idx,3),sol_stat(idx,4), sol_iter(idx,4),sol_efunc(idx,4),...
        sol_stat(idx,7), sol_iter(idx,7),sol_efunc(idx,7),sol_stat(idx,8), sol_iter(idx,8),sol_efunc(idx,8),...
        sol_stat(idx,5), sol_iter(idx,5),sol_efunc(idx,5),sol_stat(idx,6), sol_iter(idx,6),sol_efunc(idx,6),...;
        sol_stat(idx,9), sol_iter(idx,9),sol_efunc(idx,9),sol_stat(idx,10), sol_iter(idx,10),sol_efunc(idx,10));
end

end



fprintf('\n');

fprintf('                OK      MAX      FAIL   Iter       Func     F/I\n');
fprintf('Total:           %d\n', maxLen);
fprintf('       Hes_i     %d      %d       %d     %d       %d     %2.2f\n', sol_ok(1), sol_max(1), maxLen-sol_ok(1)-sol_max(1), sum(sol_iter(:,1)), sum(sol_efunc(:,1)),sum(sol_efunc(:,1))/sum(sol_iter(:,1)));
fprintf('       Hes_d     %d      %d       %d     %d       %d     %2.2f\n', sol_ok(2), sol_max(2), maxLen-sol_ok(2)-sol_max(2), sum(sol_iter(:,2)), sum(sol_efunc(:,2)),sum(sol_efunc(:,2))/sum(sol_iter(:,2)));
fprintf('       Bfgs      %d      %d       %d     %d       %d     %2.2f\n', sol_ok(3), sol_max(3), maxLen-sol_ok(3)-sol_max(3), sum(sol_iter(:,3)), sum(sol_efunc(:,3)),sum(sol_efunc(:,3))/sum(sol_iter(:,3)));
fprintf('       Sbfgs_k   %d      %d       %d     %d       %d     %2.2f\n', sol_ok(4), sol_max(4), maxLen-sol_ok(4)-sol_max(4), sum(sol_iter(:,4)), sum(sol_efunc(:,4)),sum(sol_efunc(:,4))/sum(sol_iter(:,4)));
if num_result == 10  
	fprintf('       Sbfgs_li  %d      %d       %d     %d       %d     %2.2f\n', sol_ok(7), sol_max(7), maxLen-sol_ok(7)-sol_max(7), sum(sol_iter(:,7)), sum(sol_efunc(:,7)),sum(sol_efunc(:,7))/sum(sol_iter(:,7)));
	fprintf('       Sbfgs_ld  %d      %d       %d     %d       %d     %2.2f\n', sol_ok(8), sol_max(8), maxLen-sol_ok(8)-sol_max(8), sum(sol_iter(:,8)), sum(sol_efunc(:,8)),sum(sol_efunc(:,8))/sum(sol_iter(:,8)));
end
fprintf('       Sbfgs_i   %d      %d       %d     %d       %d     %2.2f\n', sol_ok(5), sol_max(5), maxLen-sol_ok(5)-sol_max(5), sum(sol_iter(:,5)), sum(sol_efunc(:,5)),sum(sol_efunc(:,5))/sum(sol_iter(:,5)));
fprintf('       Sbfgs_d   %d      %d       %d     %d       %d     %2.2f\n', sol_ok(6), sol_max(6), maxLen-sol_ok(6)-sol_max(6), sum(sol_iter(:,6)), sum(sol_efunc(:,6)),sum(sol_efunc(:,6))/sum(sol_iter(:,6)));
if num_result == 10  
	fprintf('       Sbfgs_li  %d      %d       %d     %d       %d     %2.2f\n', sol_ok(9), sol_max(9), maxLen-sol_ok(9)-sol_max(9), sum(sol_iter(:,9)), sum(sol_efunc(:,9)),sum(sol_efunc(:,9))/sum(sol_iter(:,9)));
	fprintf('       Sbfgs_ld  %d      %d       %d     %d       %d     %2.2f\n', sol_ok(10), sol_max(10), maxLen-sol_ok(10)-sol_max(10), sum(sol_iter(:,10)), sum(sol_efunc(:,10)),sum(sol_efunc(:,10))/sum(sol_iter(:,10)));
end


for ii=1:num_result; ki(ii)=sum(sol_iter(sol_stat(:,ii)==1,ii)); end
for ii=1:num_result; kf(ii)=sum(sol_efunc(sol_stat(:,ii)==1,ii)); end

lessi = zeros(num_result,1);
lessf = zeros(num_result,1);
totalp = zeros(num_result,1);
for nSet = 4:num_result 
    for idx = 1:maxLen
		if sol_stat(idx,nSet)==1 && ( sol_stat(idx,3)==1)
			totalp(nSet) = totalp(nSet)+1;
	   		if ( sol_iter(idx,nSet) <= sol_iter(idx,3) )
	   	  		lessi(nSet) = lessi(nSet)+1;
			end
	   		if ( sol_efunc(idx,nSet) <= sol_efunc(idx,3) )
	   	  		lessf(nSet) = lessf(nSet)+1;
			end
	   	end
	end
end


fprintf('\n only for OK prob: \n');
fprintf('                OK      Iter       Func     F/I    LessI   LessF\n');
fprintf('       Hes_i      %d      %d       %d     %2.2f\n', sol_ok(1), ki(1), kf(1),kf(1)/ki(1));
fprintf('       Hes_d      %d      %d       %d     %2.2f\n', sol_ok(2), ki(2), kf(2),kf(2)/ki(2));
fprintf('       Bfgs       %d      %d       %d     %2.2f\n', sol_ok(3), ki(3), kf(3),kf(3)/ki(3));
fprintf('       Sbfgs_k    %d      %d       %d     %2.2f      %d/%d      %d/%d\n', sol_ok(4), ki(4), kf(4),kf(4)/ki(4),lessi(4), totalp(4),lessf(4), totalp(4));
if num_result == 10
fprintf('       Sbfgs_kli  %d      %d       %d     %2.2f      %d/%d      %d/%d\n', sol_ok(7), ki(7), kf(7),kf(7)/ki(7),lessi(7), totalp(7),lessf(7), totalp(7));
fprintf('       Sbfgs_kld  %d      %d       %d     %2.2f      %d/%d      %d/%d\n', sol_ok(8), ki(8), kf(8),kf(8)/ki(8),lessi(8), totalp(8),lessf(8), totalp(8));
end
fprintf('       Sbfgs_i    %d      %d       %d     %2.2f      %d/%d      %d/%d\n', sol_ok(5), ki(5), kf(5),kf(5)/ki(5),lessi(5), totalp(5),lessf(5), totalp(5));
fprintf('       Sbfgs_d    %d      %d       %d     %2.2f      %d/%d      %d/%d\n', sol_ok(6), ki(6), kf(6),kf(6)/ki(6),lessi(6), totalp(6),lessf(6), totalp(6));
if num_result == 10
fprintf('       Sbfgs_li   %d      %d       %d     %2.2f      %d/%d      %d/%d\n', sol_ok(9), ki(9), kf(9),kf(9)/ki(9),lessi(9), totalp(9),lessf(9), totalp(9));
fprintf('       Sbfgs_ld   %d      %d       %d     %2.2f      %d/%d      %d/%d\n', sol_ok(10), ki(10), kf(10),kf(10)/ki(10),lessi(10), totalp(10),lessf(10), totalp(10));
end


if num_result == 10
  num_result = 6;
  sol_iter1  = sol_iter(:,[1,2,3,4,7,8]);
  sol_efunc1 = sol_efunc(:,[1,2,3,4,7,8]);
  sol_stat1  = sol_stat(:,[1,2,3,4,7,8]);
%  do_prof_NO_INF(maxLen,num_result,sol_iter1,sol_efunc1,sol_stat1,1);
  num_result = 7;
  sol_iter1  = sol_iter(:,[1,2,3,5,6,9,10]);
  sol_efunc1 = sol_efunc(:,[1,2,3,5,6,9,10]);
  sol_stat1  = sol_stat(:,[1,2,3,5,6,9,10]);
  do_prof_NO_INF(maxLen,num_result,sol_iter1,sol_efunc1,sol_stat1,1);
else
  do_prof_NO_INF(maxLen,num_result,sol_iter,sol_efunc,sol_stat,0);
end








