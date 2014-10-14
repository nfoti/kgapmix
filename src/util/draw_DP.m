function [x] = draw_DP(C)

prob = (C - log(C)) / (C - log(C) + exp(- exp(C)));

check = 0;
while ( check == 0 )
    if ( rand < prob )
        x = C * (exp(C) / C)^rand;
        check = rand < exp(- (x - C));
    else
        x = exp(C) - log(rand);
        check = rand < exp(C) / x;
    end
end