function [D,wavespecs]=wavelet_compress(D,wavespecs,alpha)

    if nargin<3
        alpha=wavespecs.compress;
    end;
    if wavespecs.compress==1
        wavespecs.Kj_all=wavespecs.Kj;
    end;
   
    n=size(D,1);
    En=single(D);
    total_Energy=zeros(1,size(D,1));
    for i=1:n
        [Csort,ix]=sort(-abs(D(i,:)));
        total_Energy(i)=sum(Csort.^2);
        energy=cumsum(Csort.^2)/total_Energy(i);
        En(i,ix)=energy;
    end;
    %'Done Computing Relative Energy ',

    %% Now find how many coefficients are in the alpha set for at least 1,
    %%  2, 3, ..., n curves (alpha set is set of coefficients needed to 
    %%                       retain alpha proportion of total energy for 
    %%                       the function)
    Ncoeffs=zeros(size(D,1),1);
    N_alpha=sum(En<alpha);   % compute for each coefficient k, for how many
                             % of the n functions is that coefficient in
                             % the alpha set?
    %for i=1:n;    
    %    Ncoeffs(i)=sum(temp>(i-1));
    %end;
    %'Done Computing # coefficients ', toc

    %%% Now compute minimum total energy preserved if keeping the set of
    %%%  coefficients in the alpha set for at least 1,..., n curves.
    i=0;
    delta=1;
    while delta>alpha
        keep=(N_alpha>i); %%% get indicator functions of coefficients 
                            %%% in the alpha set for at least i functions
        delta=min(sum(D(:,keep)'.^2)./total_Energy);
        i=i+1;
    end;
    keep=(N_alpha>(i-2));
    temp=sum(D(:,keep)'.^2)./total_Energy;
    wavespecs.alpha_min=min(temp);
    wavespecs.alpha_mean=mean(temp);
    wavespecs.K=sum(keep);
    wavespecs.keep=keep;
    D=D(:,keep==1);
    [wavespecs.Kj,wavespecs.Kj_all]=get_Kj_compress(wavespecs);
    wavespecs.J_all=length(wavespecs.Kj_all);
    wavespecs.J=length(wavespecs.Kj);
    wavespecs.K=sum(wavespecs.Kj);
    wavespecs.compress=alpha;
end