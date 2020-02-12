%Following Algorithm 1 in section 2.4 of arXiv:1509.08391


%precision for the termination condition for the determination of the
%solution of <<x-x0-0.5delta>>=0
epsilon = 0.00001;


Emax = 55; %Maximum valueof x we want to sample
delta = 0.5; %Size of the interval
Njacknife = 20; %Number of Jacknife blocks to estimate the stochastic error of a
Naverage = 100000;%100000; %Number of random numbers we sample to calculate the reweighted expectation value <<x-x0-0.5*delta>> 




for Eint = 1:(Emax/delta) %loop through the intervals from zero to Emax with size delta

    
    %setting the lower end of the interval
    x0(Eint) = (Eint-1)*delta;
    
    a(Eint) = 0.0; %initial guess for a
    
    for jcount = 1:Njacknife %loop for Jacknife error analysis
        
        
        %initial values for a_i
        %initial value of a_i is set to 0 for the first interval, all the following are set to the value of
        %a of the previous interval
        
        if(Eint == 1)
            
            a_i(jcount) = 0.0;
        
            a_i_new = -2*epsilon;
        else
            
            a_i(jcount)=a(Eint-1) - 2*epsilon;
            
            a_i_new = a(Eint-1);
        end
        
        
        RobMarchcount = 0; %Setting the counter for the number Robson March iterations
        
        
        
        Reweightexpect = 10*epsilon;
        
        while abs(a_i_new-a_i(jcount)) > epsilon %&& RobMarchcount < 50 %termination condition for the Robson March iteration
        
            a_i(jcount) = a_i_new;
            
            %Finding right interval for uniform variable
            y1=0.5*(erf(8*a_i(jcount)+x0(Eint)/16)+1);
            y2=0.5*(erf(8*a_i(jcount)+(x0(Eint)+delta)/16)+1);
            
            
         
            %sample uniform variable in the right interval
            y = (y2-y1).*rand(1,Naverage) + y1;
                
            %sample x from the distribution
            %exp(-x^2/16^2)*exp(-ax)using en.wikipedia.org/wiki/Inverse_transform_sampling
            x = 16*erfinv(2*y-1)-128*a_i(jcount);
                
                
            %calculate reweighted expectation value <<x-x0-0.5delta>>
            Reweightexpect = sum(x)/Naverage- x0(Eint) - 0.5*delta;
            
            %Robson Monroe step
            a_i_new = a_i(jcount) + 12/(delta^2*(RobMarchcount+1))*Reweightexpect;
            %a_i_new = a_i(jcount) + Reweightexpect/var;
            RobMarchcount = RobMarchcount + 1;
            
        end
        
        
        
        a_i(jcount) = a_i_new;
        %calculate a by averaging all a_i(jcount)
        a(Eint) = a(Eint) + a_i_new/Njacknife;
        
    end
    
    s2(Eint)=sum((a_i-a(Eint)).^2)/(Njacknife-1);
    
    
    
    for i=1:Njacknife
        
        a_i_del = a_i;
        a_i_del(i) = [];
        A(Eint,i) = 1/(Njacknife-1)*sum(a_i_del);
        
    end
    
    Astar = sum(A(Eint,:))/Njacknife;
    
    var(Eint)= (Njacknife-1)*sum((Astar-A(Eint,:)).^2);
    
    errorjack(Eint) = sqrt(var(Eint))/sqrt(Njacknife);
    
    for i=1:Njacknife
        
        T = sum(A(:,i))- a(1)/2 - a(Eint)/2;
    
        rhojack(i) = exp(delta*T);
    end
    
    rho(Eint) = 1/Njacknife*sum(rhojack);
    errorrho(Eint) = sqrt((Njacknife-1)/Njacknife*sum((rho(Eint)-rhojack).^2));
    
    errordeltaa(Eint) = delta/2*rho(Eint)*sqrt(s2(Eint));
    
    
    Reweightexpect
    a(Eint)
    Eint*delta
    
end
    




% unused histogram
% x0 = 0;    
% %a = 0.0;
% 
% y1=0.5*(erf(8*a+x0/16)+1);
% y2=0.5*(erf(8*a+(x0+delta)/16)+1);
% 
% 
% for m = 1 : 1000
%     
%     r = (y2-y1).*rand() + y1;
%     
%     q = 16*erfinv(2*r-1)-128*a
%     
% end
% histogram(Q)
