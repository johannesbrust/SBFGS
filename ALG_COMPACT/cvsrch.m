      function [x,f,g,alp_trial,info,nfev,stx,sty] ...
       = cvsrch(fcn,n,x,f,g,p,alp_trial,c1,c2,xtol, ...
                 stpmin,stpmax,maxfev,prob,par,func)
%   Translation of minpack subroutine cvsrch
%   Dianne O'Leary   July 1991
%     **********
%
%     Subroutine cvsrch
%
%     The purpose of cvsrch is to find a step which satisfies 
%     a sufficient decrease condition and a curvature condition.
%     The user must provide a subroutine which calculates the
%     function and the gradient.
%
%     At each stage the subroutine updates an interval of
%     uncertainty with endpoints stx and sty. The interval of
%     uncertainty is initially chosen so that it contains a 
%     minimizer of the modified function
%
%          f(x+alp_trial*p) - f(x) - c1*alp_trial*(gradf(x)'p).
%
%     If a step is obtained for which the modified function 
%     has a nonpositive function value and nonnegative derivative, 
%     then the interval of uncertainty is chosen so that it 
%     contains a minimizer of f(x+alp_trial*p).
%
%     The algorithm is designed to find a step which satisfies 
%     the sufficient decrease condition 
%
%           f(x+alp_trial*p) <= f(x) + c1*alp_trial*(gradf(x)'*p),
%
%     and the curvature condition
%
%           abs(gradf(x+alp_trial*p)'*p)) <= c2*abs(gradf(x)'*p).
%
%     If c1 is less than c2 and if, for example, the function
%     is bounded below, then there is always a step which satisfies
%     both conditions. If no step can be found which satisfies both
%     conditions, then the algorithm usually stops when rounding
%     errors prevent further progress. In this case alp_trial only 
%     satisfies the sufficient decrease condition.
%
%     The subroutine statement is
%
%        subroutine cvsrch(fcn,n,x,f,g,p,alp_trial,c1,c2,xtol,
%                          stpmin,stpmax,maxfev,info,nfev,wa)
%     where
%
%	fcn is the name of the user-supplied subroutine which
%         calculates the function and the gradient.  fcn must 
%      	  be declared in an external statement in the user 
%         calling program, and should be written as follows.
%
%         function [f,g] = fcn(x,prob,par,func) (Matlab)     (10/2010 change in documentation)
%	  (derived from Fortran subroutine fcn(n,x,f,g) )
%         integer n
%         f
%         x(n),g(n)
%	  ----------
%         Calculate the function at x and
%         return this value in the variable f.
%         Calculate the gradient at x and
%         return this vector in g.
%	  ----------
%	  return
%	  end
%
%       n is a positive integer input variable set to the number
%	  of variables.
%
%	x is an array of length n. On input it must contain the
%	  base point for the line search. On output it contains 
%         x + alp_trial*p.
%
%	f is a variable. On input it must contain the value of f
%         at x. On output it contains the value of f at x + alp_trial*p.
%
%	g is an array of length n. On input it must contain the
%         gradient of f at x. On output it contains the gradient
%         of f at x + alp_trial*p.
%
%	p is an input array of length n which specifies the
%         search direction.
%
%	alp_trial is a nonnegative variable. On input alp_trial contains an
%         initial estimate of a satisfactory step. On output
%         alp_trial contains the final estimate.
%
%   c1 and c2 are nonnegative input variables. Termination
%         occurs when the sufficient decrease condition and the
%         directional derivative condition are satisfied.
%
%	xtol is a nonnegative input variable. Termination occurs
%         when the relative width of the interval of uncertainty 
%	  is at most xtol.
%
%	stpmin and stpmax are nonnegative input variables which 
%	  specify lower and upper bounds for the step.
%
%	maxfev is a positive integer input variable. Termination
%         occurs when the number of calls to fcn is at least
%         maxfev by the end of an iteration.
%
%	info is an integer output variable set as follows:
%	  
%	  info = 0  Improper input parameters.
%
%	  info = 1  The sufficient decrease condition and the
%                   directional derivative condition hold.
%
%	  info = 2  Relative width of the interval of uncertainty
%		    is at most xtol.
%
%	  info = 3  Number of calls to fcn has reached maxfev.
%
%	  info = 4  The step is at the lower bound stpmin.
%
%	  info = 5  The step is at the upper bound stpmax and satisfies Wolfe
%	  conditions
%
%	  info = 6  Rounding errors prevent further progress.
%                   There may not be a step which satisfies the
%                   sufficient decrease and curvature conditions.
%                   Tolerances may be too small.
%
%	  info = 7  The step is at the upper bound stpmax and does not satisfy
%   Wolfe conditions
%
%       nfev is an integer output variable set to the number of
%         calls to fcn.
%
%	wa is a work array of length n.
%
%     Subprograms called
%
%	user-supplied......fcn
%
%	MINPACK-supplied...cstep
%
%	FORTRAN-supplied...abs,max,min
%	  
%     Argonne National Laboratory. MINPACK Project. June 1983
%     Jorge J. More', David J. Thuente
%
%     **********
      p5 = .5;
      p66 = .66;
      xtrapf = 4;
      info = 0;
      infoc = 1;
      nfev = 0;
      
      stx = 0.0;
      sty = 0.0;      
%
%     Check the input parameters for errors.
%
      if (n <= 0 | alp_trial <= 0.0 | c1 < 0.0 |  ...
          c2 < 0.0 | xtol < 0.0 | stpmin < 0.0  ...
          | stpmax < stpmin | maxfev <= 0) 
         return
      end
%
%     Compute the initial gradient in the search direction
%     and check that p is a descent direction.
%
      dginit = g'*p;
      if (dginit >= 0.0) 
          return
      end
%
%     Initialize local variables.
%
      brackt = 0;
      stage1 = 1;
      finit = f;
      dgtest = c1*dginit;
      width = stpmax - stpmin;
      width1 = 2*width;
      wa = x;
%
%     The variables stx, fx, dgx contain the values of the step, 
%     function, and directional derivative at the best step.
%     The variables sty, fy, dgy contain the value of the step,
%     function, and derivative at the other endpoint of
%     the interval of uncertainty.
%     The variables alp_trial, f, dg contain the values of the step,
%     function, and derivative at the current step.
%

      
      fx = finit;
      dgx = dginit;
      
      fy = finit;
      dgy = dginit;
%
%     Start of iteration.
%
    maxiter = 1000;
    num_iter = 0;
    while (num_iter < maxiter)   
%
%        Set the minimum and maximum steps to correspond
%        to the present interval of uncertainty.
%
         if (brackt) 
            stmin = min(stx,sty);
            stmax = max(stx,sty);
         else
            stmin = stx;
            stmax = alp_trial + xtrapf*(alp_trial - stx);
         end 
%
%        Force the step to be within the bounds stpmax and stpmin.
%
         alp_trial = max(alp_trial,stpmin);
         alp_trial = min(alp_trial,stpmax);
%
%        If an unusual termination is to occur then let 
%        alp_trial be the lowest point obtained so far.
%
         if ((brackt & (alp_trial <= stmin | alp_trial >= stmax)) ...
            | nfev >= maxfev-1 | infoc == 0 ...
            | (brackt & stmax-stmin <= xtol*stmax)) 
            alp_trial = stx;
         end
%
%        Evaluate the function and gradient at alp_trial
%        and compute the directional derivative.
%
         x = wa + alp_trial * p;
         [f,g, gradu,Hkplus] = feval(fcn,x,prob,par,func);
         nfev = nfev + 1;
         dg = g' * p;
         ftest1 = finit + alp_trial*dgtest;
%
%        Test for convergence.
%
         if ((brackt & (alp_trial <= stmin | alp_trial >= stmax)) | infoc == 0) 
                  info = 6;
         end
         if (alp_trial == stpmax & f <= ftest1 & dg <= dgtest) 
             if abs(dg) <= c2*(-dginit)
                  % satisfies Wolfe conditions
                  info = 5;
             else
                 info = 7;
             end
         end
         if (alp_trial == stpmin & (f > ftest1 | dg >= dgtest)) 
                  info = 4;
         end
         if (nfev >= maxfev) 
                  info = 3;
         end
         if (brackt & stmax-stmin <= xtol*stmax) 
                  info = 2;
         end
         if (f <= ftest1 & abs(dg) <= c2*(-dginit)) 
              info = 1;
         end
         
         % do modified LS only for sbfgs-m with full rank update
         if strcmp(par.alg,'s-bfgs')==1 && par.sbfgsAlg==1 && strcmp(par.addUnknown,'addNonLinObj_Proj')==0 && ((info==5) || (info ==1 ))
            % Wolfe conditions satisfied
            %test curvature condition: bar_gamma > 0
            
            % 04/26/19, J.B. count function evaluations
            nfev = nfev + 1;
            
            [f0,g0, gradu0,~] = feval(fcn,wa,prob,par,func);
            sk = x - wa;
            yk = gradu-gradu0;
            
            if strcmp(par.addUnknown,'addNonLinObj_Proj')
                Pmat = prob.Pmat;
                sk = Pmat*sk;
                Hkplus = Pmat*Hkplus*Pmat';
            end
            
            bk = yk + Hkplus*sk;
            gamma  = sk'*bk; 
            if gamma<=0
                info = 0;
            end         
         end         
         
         %
%        Check for termination.
%
         if (info ~= 0) 
                  return
         end
%
%        In the first stage we seek a step for which the modified
%        function has a nonpositive value and nonnegative derivative.
%
         if (stage1 & f <= ftest1 & dg >= min(c1,c2)*dginit) 
                stage1 = 0;
         end
%
%        A modified function is used to predict the step only if
%        we have not obtained a step for which the modified
%        function has a nonpositive function value and nonnegative 
%        derivative, and if a lower function value has been  
%        obtained but the decrease is not sufficient.
%
         if (stage1 & f <= fx & f > ftest1) 
%
%           Define the modified function and derivative values.
%
            fm = f - alp_trial*dgtest;
            fxm = fx - stx*dgtest;
            fym = fy - sty*dgtest;
            dgm = dg - dgtest;
            dgxm = dgx - dgtest;
            dgym = dgy - dgtest;
% 
%           Call cstep to update the interval of uncertainty 
%           and to compute the new step.
%
            [stx,fxm,dgxm,sty,fym,dgym,alp_trial,fm,dgm,brackt,infoc] ...
             = cstep(stx,fxm,dgxm,sty,fym,dgym,alp_trial,fm,dgm, ...
                     brackt,stmin,stmax);
%
%           Reset the function and gradient values for f.
%
            fx = fxm + stx*dgtest;
            fy = fym + sty*dgtest;
            dgx = dgxm + dgtest;
            dgy = dgym + dgtest;
         else
% 
%           Call cstep to update the interval of uncertainty 
%           and to compute the new step.
%
            [stx,fx,dgx,sty,fy,dgy,alp_trial,f,dg,brackt,infoc] ...
             = cstep(stx,fx,dgx,sty,fy,dgy,alp_trial,f,dg, ...
                     brackt,stmin,stmax);
         end
%
%        Force a sufficient decrease in the size of the
%        interval of uncertainty.
%
         if (brackt) 
            if (abs(sty-stx) >= p66*width1) 
              alp_trial = stx + p5*(sty - stx);
            end
            width1 = width;
            width = abs(sty-stx);
         end
         
         num_iter = num_iter+1;
%
%        End of iteration.
%
     end
%
%     Last card of subroutine cvsrch.
%


