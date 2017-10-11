
        
function[sigmak] = mise_a_jour_sigma(lm_pq,sigmak,sigmabar)


if (sigmak < max(lm_pq)+sigmabar)
			sigmak = max(1.5*sigmak,max(lm_pq)+sigmabar);
		else
			if (sigmak > 1.1*(max(lm_pq)+sigmabar))
				sigmak = (sigmak+max(lm_pq)+sigmabar)/2;
            else sigmak = sigmak ;
            end 
        end


end
