function [e,ce,g,ae,hl,indics] = chs (indic,xy,lmde) 
    global LB
    global a
    global b
    n = length(xy)/2 ;
    nb = n + 1 ;
    
    e = 0;
    ce = sparse(nb,1);
    g = sparse(2*n,1);
    ae = sparse(nb,2*n);
    hl = sparse(2*n,2*n);
    
    switch indic 
        case 1
            dessine_chaine(xy,lmde);
            
        case 2
            e = calcul_energie(xy,lmde);
            ce =calcul_contrainte_egal(xy,lmde);
        case 4
            e = calcul_energie(xy,lmde);
            ce =calcul_contrainte_egal(xy,lmde);
            g = calcul_gradient(xy,lmde);
            ae = calcul_jacobien(xy,lmde);

        case 5
            hl = calculderiv_num_ieme_x_hessien(xy,lmde);
            
        otherwise % problems
            indics = 1;
    end
    
    indics = 0;
    
    
    
end
    
%functions to get every parameter

function [] = dessine_chaine(xy,lmde)

  global LB
   % global (a,b)
    global a
    global b
n = length(xy)/2 ;
x = [0; xy( 1:n) ; a];
y = [0; xy(n+1: end) ;  b] ;
plot(x,y);
%add plot 
end



function [e] = calcul_energie(xy,lmde)
    global LB
   % global (a,b)
    global a
    global b
    
    n = length(xy)/2 ;
    y = [ 0 ; xy(n+1 : end); b];
    milieux = (y(2:end)+y(1:end-1))/2 ;
    y = [];
    e = LB'*milieux ;

end

function[ce] = calcul_contrainte_egal(xy,lmde)
    global LB
   % global (a,b)
    global a
    global b
    n = size(xy,1)/2 ;
    nb = n + 1 ;
    
  ce(1)=xy(1)^2+xy(n+1)^2-LB(1)^2 ;
  ce(2:n)=xy(2:n).*xy(2:n)+xy(1:n-1).*xy(1:n-1)+xy(n+2:2*n).*xy(n+2:2*n)+xy(n+1:2*n-1).*xy(n+1:2*n-1)-2*xy(2:n).*xy(1:n-1)-2*xy(n+2:2*n).*xy(n+1:2*n-1)-LB(2:n).*LB(2:n);
  ce(nb)=(a-xy(n))^2+(b-xy(2*n))^2-LB(nb)^2 ;
  
  ce;
  
  ce = ce';
    
   % n = length(xy)/2;
   % y = [0 ; xy(n+1 : end); b];
   % x = [0; xy( 1:n) ; a];
    
   % vecteur2 = ((y(2:end)-y(1:end-1)));
   % vecteur1 = ((x(2:end)-x(1:end-1)));
    
   % ce = vecteur1.^2+vecteur2.^2- LB.^2
    

end


function [g] = calcul_gradient(xy,lmde)

   global LB
   % global (a,b)
   global a
   global b
   n=length(xy)/2 ;
   g = [zeros(n,1); (LB(1:end-1)+LB(2:end))./2];
   
end 




function [DF , pas] = gradient_numerique(f,x,lmde)

n = length(x);
DF = zeros(n,1);
pas = zeros(n,1);
for i=1:n
    
    [DF(i,1), pas(i)]=deriv_num_ieme_x(f,i,x,lmde) ;
    
end    

end

function [y, pas] = deriv_num_ieme_x(f,i,x,lmde)

tau = 1 ;
pas = sqrt(eps)*max(tau,abs(x(i)));
ei = zeros(size(x));
ei(i) = 1;

y = ((calcul_energie(x+pas*ei,lmde))-(calcul_energie(x-pas*ei,lmde)))/(2*pas) ;

end

function [ae] = calcul_jacobien(xy,lmde)
   global LB
   % global (a,b)
   global a
   global b
   n=length(xy)/2 ;
  
   
%Jacobian ce

   %Step 1 : son variables x 
   %work only on consecutif Der !!!
   %V1 dfi/xi
   v1(1)=2*(xy(1)-0);
   v1(2:n)=2*(xy(2:n)-xy(1:n-1));
   %V2 dfi/x(i+1)
   v2(1:n-1)=-2*(xy(2:n)-xy(1:n-1));
   v2 = [v2,0] ;% add 0
   aex=spdiags([v1',v2'],[0,-1],n,n);
   %add last element
   aex(n+1,n)=-2*(a-xy(n));
   
   %Etape2 : on  y from n+1 to 2n 
   
   
   w1(1)=2*(xy(n+1)-0);
   w1(2:n)=2*(xy(n+2:2*n)-xy(n+1:2*n-1));
   %w2 dfi/y(i+1)
   w2(1:n-1)=-2*(xy(n+2:2*n)-xy(n+1:2*n-1));
   w2 = [w2,0]; % rajouter 0
   aey=spdiags([w1',w2'],[0,-1],n,n);
   %rajouter dernier element
   aey(n+1,n)=-2*(b-xy(2*n));
   
   
    ae=[aex,aey];
end


%hl

function [hl] = calculderiv_num_ieme_x_hessien(xy,lmde)
 global LB
   % global (a,b)
   global a
   global b
   n=length(xy)/2 ;
   nb = n +1 ;
   
  hl=sparse(2*n) ;
  %Get the diagonal
  hl(1:n-1,1:n-1)=2*(diag(lmde(1:n-1)+lmde(2:n))-diag(lmde(2:n-1),1)-diag(lmde(2:n-1),-1));
  %extremity values
  hl(n,n)=2*(lmde(nb)+lmde(n));
  if n-1 >0 
      hl(n,n-1)=-2*lmde(n) ;
      hl(n-1,n)=-2*lmde(n);
  end
  %Complete the Hessian
  hl(n+1:2*n-1,n+1:2*n-1)=2*(diag(lmde(1:n-1)+lmde(2:n))-diag(lmde(2:n-1),1)-diag(lmde(2:n-1),-1)) ;
  hl(2*n,2*n)=2*(lmde(nb)+lmde(n));
  if n-1>0
      hl(2*n,2*n-1)=-2*lmde(n) ;
      hl(2*n-1,2*n)=-2*lmde(n);
  end
  

%hl
end




