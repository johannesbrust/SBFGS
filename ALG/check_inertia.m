function [Mat_new,par]=check_inertia(n,m,Mat,rhs,prob,par_in,numIter)

pri_reg  = 0;
dual_reg = 0;

par = par_in;

alg = par.checkInertia;

in_done = 0;

if (alg==1)
    % check inertia
    eK = eig(Mat);
    npos = length(find(eK > 0));
    nneg = length(find(eK < 0));
    par.numFact = par.numFact + 1;
    % if inertia is correct, return
    if (npos==n && nneg==m)
        par.done_correction = 0;
        Mat_new = Mat;
        return;
    end
elseif alg==2
    % check descent
    grad_0		= prob.obj_grad;
    p = Mat\rhs;
    if grad_0'*p <0
        par.done_correction = 0;
        Mat_new = Mat;
        return;
    end    
end

% start to add regularization
initReg = par.init_reg; % initial value of primal reg
reg_old = par.last_reg;	% previous value of primal reg

% linear solver params   
       		trial 	= 1;        % inertia trial counter
   		step_done 	= 0;        % step flag

use_specialInitVal = 0;
while step_done == 0  && trial <= 50
 	% add regularization 
    if trial == 1
       if reg_old == 0
           	pri_reg = initReg;
       else
           	pri_reg = reg_old/3;
       end
       
       if numIter > 0 && strcmp(par.alg,'s-bfgs')&& ( (par.sbfgsAlg==2) && strcmp(par.addUnknown,'addNonLinObj')==1 )
           if strcmp(par.addUnknown,'addNonLinObj_Proj')
               sk = prob.sk;
               Hes = prob.Pmat*prob.Hes*prob.Pmat';
               yk = prob.yk;
               initVal4sbfgsplus = (1e-5-((yk+Hes*sk)')*sk)/norm(sk,2);
           else
               initVal4sbfgsplus = (1e-5-((prob.yk+prob.Hes*prob.sk)')*prob.sk)/norm(prob.sk,2);
           end
           if pri_reg < initVal4sbfgsplus
               pri_reg = initVal4sbfgsplus;
               use_specialInitVal = 1;
           end
       end
    else
		if reg_old == 0
       		pri_reg = min(pri_reg*100,1e+20);
       	else
       		pri_reg = min(pri_reg*8,1e+20);
       	end
    end
    
%     if (npos+nneg~=n+m)      % if system is singular
%         dual_reg = 1e-10;    
%     end
        
	diagMat = pri_reg*eye(n);
	if m>0
		diagMat = [diagMat 			zeros(n,m)
				   zeros(m,n)	  -dual_reg*eye(m)]
    end
	Mat_new =  Mat + diagMat;

    if (alg==1)
        % check inertia    
        % get number of positive and negative eigs
        eK = eig(Mat_new);
        npos = length(find(eK > 0));
        nneg = length(find(eK < 0));
        par.numFact = par.numFact + 1;
        
        % if inertia is correct, return
        if (npos==n && nneg==m)
            par.done_correction = 1;
            par.last_reg = pri_reg;
            if use_specialInitVal==1 && trial == 1
                warning('Use Special initial value for regularization and pass the inertia test in the 1st iter');
            end
            return;
        else
            trial = trial + 1;
            use_specialInitVal = 0;
        end
    elseif alg==2  
        % check descent
        p = Mat_new\rhs;
        if grad_0'*p <0
            par.done_correction = 1;
            par.last_reg = pri_reg;
            if use_specialInitVal==1 && trial == 1
                warning('Use Special initial value for regularization and pass the inertia test in the 1st iter');
            end            
            return;
        else
            trial = trial + 1;
            use_specialInitVal = 0;
        end  
    end
end

if(trial > 50)
    warning('Reg fail: primal reg becomes too large');
end

