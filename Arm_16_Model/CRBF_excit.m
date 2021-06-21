function [excitation] =CRBF_excit(time,muscle)
% written by Russell T Johnson, University of Southern California
% rtjohnso@usc.edu
% Last edited: 6/18/2021    

% calculate the muscle excitation signal for a single muscle based on
% compact radial basis functions. There are 10 functions. 
    
    nNodes = 10; 
    interval = 0.600/(nNodes-1); 
    c = -0.05:interval:0.55;
    w = interval*1.33;

    for j = 1:nNodes
        output(:,j) = muscle(j).*bumpfun(time,c(j),w);
    end

    raw_out = sum(output,2); % +constant; There's an option to add a constant to shift, but I think as long as the parameters can be negative, we don't need these values 

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
