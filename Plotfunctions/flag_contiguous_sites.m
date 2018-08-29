function flag2=flag_contiguous_sites(flag)

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
            
        

