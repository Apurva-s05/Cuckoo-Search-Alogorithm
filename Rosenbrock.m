% ==================================================================================================
%                             Rosenbrock Function-- MULTIMODAL
% ==================================================================================================

%Phase 1- INPUT PARAMETERS 
D=100; 
lb=-5*ones(1,D); % Lower bound  
ub=5*ones(1,D);  % Upper bounds
N = 20;         %population size
n = N; 
pa = 0.25;      %discovery rate of alien eggs/solution
max_iter = 10; %maximum iterations
 
%Phase 2- Defining objective fuction 
   for i=1:N,
      nest(i,:)=lb+(ub-lb).*rand(size(lb));
   end
fx = fns(nest);
beta = 3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

% PHASE 3: CUCKOO SEARCH MAIN LOOP START
for iter =1:max_iter
    [fnv,indf] = min(fx);    %fnv is function value and indf is index value
    best = nest(indf,:);
    
    for j = 1:N
        s = nest(j,:);
        X=s;
        
 %Phase 4- Levy flights by Mantegna's algorithm
    u=randn(size(s))*sigma;
    v=randn(size(s));
    step=u./abs(v).^(1/beta);
    Xnew = X+randn(size(s)).*0.01.*step.*(X-best); 
 
 
% Phase 5: check bounds
 for kk = 1:size(Xnew,2)
     if Xnew(kk)>ub(kk)
         Xnew(kk) = ub(kk);
     elseif Xnew(kk)<lb(kk)
          Xnew(kk) = lb(kk);
     end
 end 
    
% Perform greedy selection
   fnew = fns(Xnew);
   if fnew<fx(j,:)
       nest(j,:) = Xnew;
       fx(j,:) = fnew;
   end
    end
 
% find the current best
    [fmin,K1] = min(fx);
    best = nest(K1,:);
    
    
% REPLACE SOME NEST BY CONSTRUCTING NEW SOLUTION
    K = rand(size(nest))< pa;
    stepsizeK = rand*(nest(randperm(n),:)-nest(randperm(n),:));
    new_nest = nest+stepsizeK.*K;
    
% check bounds
for ii=1:size(nest,1)
    s = new_nest(ii,:);
    for kk = 1:size(s,2)
    if s(kk)>ub(kk)
       s(kk) = ub(kk);
    elseif s(kk)< lb(kk)
            s(kk) = lb(kk);
       end
   end
    new_nest(ii,:)=s;
    
% Perform greedy selection
    fnew = fns(s);
    if fnew<fx(ii,:)
        nest(ii,:)=s;
        fx(ii,:)=fnew;
    end
    
end

% Memorize the best
[optval,optind] = min(fx);
BestFx(iter)=optval;
BestX(iter,:)=nest(optind,:);
 
% Ploting the result
plot(BestFx,'Linewidth',2);
xlabel('Iteration Number');
ylabel('Fitness value')
title('Convergence Vs Iteration')
grid on
 
 
end


% Objective function
function out = fns(X) 
 out = 0;
    n = size(X, 2);
    assert(n >= 1, 'Given input X cannot be empty');
    a = 1;
    b = 100;
    for i = 1 : (n-1)
        out = out + (b * ((X(:, i+1) - (X(:, i).^2)) .^ 2)) + ((a - X(:, i)) .^ 2);
    end
end
