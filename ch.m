function [res] = ch()
%Donnees Initial 
    global LB
    global a
    global b
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Choise of optimization algorithm
% imode = 1 : Newton,
% imode =0 : Newton avec RL

    imode = 0;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   %Information on max-iteration and tolerence :
   
    maxit=50;
    tol=[10^(-6) 10^(-6)];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Information on initial chain
    % cas1 - 
    
    %LB = [0.7 0.5 0.3 0.2 0.5]'; 
    %a = 1;
    %b=-1;
    
    %xy=[0.2 0.4 0.6 0.8 -1 -1.5 -1.5 -1.3]'; 
    a = 1 ;
    b = -0.3;
    LB = [0.6 0.4 0.2 0.2 0.4]';
    xy = [0.5  0.8  1.1  1.3 ...
        -0.5 -0.6 -0.6 -0.6]';

    
   % cas  -
   
    %LB = [0.7 0.5 0.3 0.2 0.5]'; 
    %a = 1;
    %b=-1;
    %xy=[0.2 0.4 0.6 0.8 1 1.5 1.5 1.3]'; 

   %cas  - 
   
    %LB = [0.7 0.5 0.3 0.2 0.5]'; 
    %a = 1;
    %b=-1;
    %xy=[0.2 0.4 0.6 0.8 -1 -1.5 1.5 -1.3]'; 

    
   % cas - 
   
    %LB = [0.7 0.5 0.3 0.2 0.5]'; 
    %a = 1;
    %b=-1;
    %xy=[0.2 0.4 0.6 0.8 1 1.5 1.5 1.3]'; 

    
   % cas 5 - 
   
   %LB = [0.5 0.5 2.0 0.4 0.4]';
   %a = 0 ;
   %b = -1;
   %xy = [ 0.16  0.38 -0.20 -0.10 0.12  0.45 -1.60 -1.30 ]';
      
%

   % xy = [1.5 -0.5 ]';
   % LB = [1 1]';
   % a = 2 ;
   % b = 0 ;
  
    %xy = [0.2 0.4 0.6 0.8 1 1.5 1.5 1.3]';
    %xy = [0.2 0.4 0.6 0.8 1 -1.5 1.5 -1.3]';
    %xy = [0.2 0.4 0.6 0.8 -1 -1.5 1.5 -1.3]';
    
    
    %LB = [0.6 0.4 0.2 0.2 0.4]';
    %xy = [0.5  0.8  1.1  1.3 ...
    %  -0.5 -0.6 -0.6 -0.6]';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Call function chs to calculate e, ce, g, ae, h1, indics(=1 errors)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %lmde = [0 0]' ;
    lmde=[0 0 0 0 0]';
    [e, ce , g , Ae , ~ , indics] = chs(4, xy , lmde);
    [ ~ , ~ , ~ , ~  , hl , indics] = chs(5, xy , lmde) ;
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Initialisation :
    
    %lmde = [0.5077 0.4223 0.5190 0.6156 0.8774]';
    %lmde = [  1 ;1; 0.41;3;3];
    lmde = -(Ae*Ae')\(Ae*g); 
    [e, ce , g , Ae , ~ , indics] = chs(4, xy , lmde);
    [ ~ , ~ , ~ , ~  , hl , indics] = chs(5, xy , lmde) ;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Save results 
    
    Nom_fichiertext = 'Results.txt';
    fout = fopen(Nom_fichiertext,'w'); 
    
    if(imode == 0)
        fprintf(1,'%s \n\n','!!! Algorithme de Newton General !!!!');
        
    else 
        fprintf(1,'%s \n\n','!!! Algorithme de Newton !!!!');

        fprintf(fout,'%s \t \t %s \t\t\t %s \t\t %s \t\t %s \n','Nb iter','f','|Grad_Lagr|','|ce|','|Lambda|');
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply the optimizer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   fprintf(1,'%s \n %s \n %s','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%','Please see optimization results in results.txt','%%%%%%%%%%%%%%%%%%%%%%%%%%%')  
    
   sqp('chs',xy,lmde,e,ce,g,Ae,hl,tol,maxit,imode,fout);
    
    
    
end    

