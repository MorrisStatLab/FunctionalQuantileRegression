function flag2=flag_contiguous_sites(flag)

% Only regions with at least 3 consecutive locations flagged based on 
% SimBaS and minimum effect size are considered significant.
% This is done to avoid flagging of singletons which are
% likely to be false positives.
flag2=zeros(size(flag));

flag=[0 0 flag 0 0];

for i=1:size(flag2,2)
    if flag(i+2)==1 && flag(i+1)==1 && flag(i)==1
        flag2(i)=1;
    elseif flag(i+2)==1 && flag(i+1)==1 && flag(i+3)==1
        flag2(i)=1;
    elseif flag(i+2)==1 && flag(i+3)==1 && flag(i+4)==1
        flag2(i)=1;
    end
end
            
        

