function [x,lm,mode] = sqp (simul, x,lm,f,c,g,A,L, tol, maxit,imode,fout) %fout,imode)

iter = 0 ;
n = size(x,1);
m=size(c,1);
iter=0; %Iteration counter
simu=1; %Simulation counter
Vitesse = [0]; % Get the speed of convergance
%Val_f = [] ; 
norme_gradient_lagrang_1 = norm(g+A'*lm , 'inf') ;
norm_contraintes_1 = norm(c,'inf') ;
omega = 0.00001;
grad= g+A'*lm;
sigmabar=max(sqrt(eps),max(lm)/100);
sigmak=max(lm)+sigmabar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%imode == 0 : Use the Global Method (SQP and Line search)


if (imode == 0)

while 1  
         
        [LL,DD] = cholmod(L,1e-6,1e23); %  
        M = sparse(LL*diag(DD)*LL');     
        AA = [M, A' ; A, sparse(m,m)] ;        
        delta = -AA\[g;c];
        dk = delta(1:n);
        lm_pq = delta(n+1:end);
        
        % Calculate sigma_k

        [sigmak] = mise_a_jour_sigma(lm_pq,sigmak,sigmabar)   ;     
            
        % calculate alpha   
       
        [alpha,simu,pas,decroissance,estimpente] = calcule_alpha(sigmak,simul,simu,x,dk,lm,M,omega,c,f,lm_pq,fout);       
    
        %test de convergance 
        
        grad= g+A'*lm;
        %grad = 0 ;
 
        norme_gradient_lagrang = norm(grad , 'inf') ;
        norm_contraintes = norm(c,'inf') ;
        conv = (norme_gradient_lagrang < tol(1) && norm_contraintes < tol(2)) ;
  
            if(conv)
              mode = 0
              break 
            end
            
        % Update
        iter = iter +1 ;

        x = x + alpha*dk;
    	lm = lm + alpha*(lm_pq-lm);
        [f, c , g , A , ~ , indics] = chs(4, x , lm);
        [ ~ , ~ , ~ , ~  , L , indics] = chs(5, x , lm) ;
        %Iteration test
         
         if iter==maxit+1
                 mode = 2 
                break  ;
         end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Save Results
    


    fprintf(fout,'%s  %s \t %s \t\t\t %s \t\t %s \t\t %s  \t\t %s \n ','N iter','Nb Simul','f','|Grad_Lagr|','|ce|', 'alpha','|Lambda_Pq|');
    fprintf(fout,'%d \t \t %d \t\t %e \t %e \t %e \t %e \t %e \n\n',[iter,simu,f,norm(g+A'*lm,'inf'),norm(c,'inf'), alpha,norm(lm,'inf') ]) ;
    fprintf(fout,'%s\n\n','Recherch Lin Armijo');
    fprintf(fout,'%s \t \n','--- Pas ---');
    fprintf(fout,' %i\t  \n',pas);
    fprintf(fout,'%s \t \n','-Decroissance-');
    fprintf(fout,' %i\t  \n',decroissance);
    fprintf(fout,'%s \t \n','--Estim-pente--');
    fprintf(fout,' %i\t  \n',estimpente);
    fprintf(fout,'%s\n\n','---------------------------------------');
    
  %%%%%%%%%% Speed
   grad=g+A'*lm;
   Vitesse=[Vitesse,log((norm(grad,'inf')+norm(c,'inf'))/(norme_gradient_lagrang_1+norm_contraintes_1))];
       
end

   %plot
    feval(simul,1,x,lm);
    %plot_vitesse
    %plot(Vitesse);    
x
if (iter <= maxit+1)
%Nature of solution : 
    r=eig(L);
	if mode ==0
        
	if max(r)<0
   	disp('c''est un maximum local')
	elseif min(r)>0
   	disp('c''est un minimum local')
	else
   	disp('un point selle')
    end 
    end

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%% imode == 1 - Methode Newton

if (imode == 1)
    
while 1 % 
  
    AA = [L,A'; A , sparse(zeros(size(A,1),size(A',2)))];
    V = -([g;c]) ; 
    deltaZ = AA\V ;    
    %Update 

    x = x + deltaZ(1:n);
    lm = deltaZ(n+1:end);
    [f, c , g , A , ~ , indics] = chs(4, x , lm);
    [ ~ , ~ , ~ , ~  , L , indics] = chs(5, x , lm) ;
    
    if(indics == 1)
        mode = 1
        break;
    end

    iter = iter +1;
    
    if iter==maxit+1
             mode = 2 
             break  ;
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%Save Results      
    fprintf(fout,'%d  \t \t %e \t \t %e \t %e  \t \t %e \n\n\n',[iter,f,norm(g+A'*lm,'inf'),norm(c,'inf'),norm(lm,'inf')]);
  
     %test de convergance 
    grad= g+A'*lm;
    %grad = 0 ;
        
    norme_gradient_lagrang = norm(grad , 'inf') ;
    norm_contraintes = norm(c,'inf') ;
    
    conv = (norme_gradient_lagrang < tol(1) && norm_contraintes < tol(2)) ;
    
    if(conv)
        mode = 0
        break 
    end
%%%%%%%%%%%Speed

    grad=g+A'*lm;
    Vitesse=[Vitesse,log((norm(grad,'inf')+norm(c,'inf'))/(norme_gradient_lagrang_1+norm_contraintes_1))];
        
end

x 
feval(simul,1,x,lm);
%Vitesse
%plot(Vitesse);

%Nature of the solution : 

if (iter <= maxit+1)
%Nature de la solution : 
    r=eig(L);
	if mode ==0
        
	if max(r)<0
   	disp('c''est un maximum local')
	elseif min(r)>0
   	disp('c''est un minimum local')
	else
   	disp('un point selle')
    end 
    end

end
    

end

end


 
   