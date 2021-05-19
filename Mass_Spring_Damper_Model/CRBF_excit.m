function [excitation] =CRBF_excit(time,muscle)
    % calculate the muscle excitation signal for a single muscle based on
    % compact radial basis functions. There are 8 functions. 

    nNodes = 10; 
    interval = 0.600/(nNodes-1); 
    c = -0.05:interval:0.55;
    w = interval*1.33;

    for j = 1:nNodes
        output(:,j) = muscle(j).*bumpfun(time,c(j),w);
    end

    raw_out = sum(output,2); % 

    for i = 1:length(raw_out)
        if raw_out(i)>0 
            excitation(i) = 1./(1+exp(-raw_out(i)));
        else
            excitation(i) = exp(raw_out(i))./(1+exp(raw_out(i)));
        end
    end

end

function output = bumpfun(time,c,s)
    for i = 1: length(time)
        if(abs((time(i)-c)/s) > 1.0)
            output(i) = 0.0;
        else
            output(i) = (exp(1.0-1.0/(1.0-((time(i)-c)/s).^2)));
        end
    end
    
end
