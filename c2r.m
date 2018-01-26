%% Function to turn a column vector into a row vector
function[row] = c2r(column)
    
    [m,n] = size(column);
    standin = ones(n);
    
    i = 1;
    
    while i < m;
        
        standin(i) = column(m-i);
        i = i+1;
        
    end
    
row = standin;
end