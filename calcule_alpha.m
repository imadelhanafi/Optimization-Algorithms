      

function [alpha,simu,pas,decroissance,estimpente] = calcule_alpha(sigmak,simul,simu,x,dk,lm,M,omega,c,f,lm_pq,fout)

		alpha=1;
        theta_init = f+sigmak*norm(c,1);
        delta0 = -dk'*M*dk + lm_pq'*c - sigmak*norm(c,1);
        
        %(-dk'*LL*diag(DD)*LL'*dk+lm'*c-sigmak*norm(c,1))
		[ff,cc,~,~,~,indics]=feval(simul,4,x+ alpha*dk,lm);
         %iter = iter +1 ;
         simu=simu+1;
       
        theta = ff+sigmak*norm(cc,1) ;
       
       
      % conver = theta <= theta0 + omega*alpha*delta0 % cri de convergance
    pas=[1]; %pas
    estimpente = [theta-theta_init]; % la pente
    decroissance=[(theta-theta_init)/alpha]; %decroissance
    
         
		%while (theta - (theta0 + omega*alpha*delta0) > 0)
         while    theta - (theta_init + omega*alpha*delta0) > 0
             
            if(indics == 1)
            mode = 1
            break;
            end
        
			alpha = alpha/2;
			[ff,cc,~,~,~,indics]=feval(simul,4,x+alpha*dk,lm);
            simu=simu+1;
            theta = ff+sigmak*norm(cc,1);
            
            pas=[pas alpha];
            decroissance=[decroissance theta-theta_init];
            estimpente=[estimpente (theta-theta_init)/alpha];
            

         end

        
end