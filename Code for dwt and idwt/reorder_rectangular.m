function [reorder,Kj]=reorder_rectangular(Kj1,Kj2)

%%%% [reorder]=reorder_rectangular(Kj1,Kj2);
%%%% Reorder wavelet coefficients from a rectangular 2d wavelet transform
%%%% with Kj1 a vector of length J1+1 containing the number of wavelet
%%%% coefficients per level for 1d wavelet transform applied to the rows,
%%%% and Kj2 a vector of length J2+1 containing the number of wavelet
%%%% coefficients per level for 1d wavelet transform applied to the
%%%% columns. Note: be careful to remove the last element of Kj1 and Kj2,
%%%% which are T1 and T2, respectively.
%%%% where K= K1*K2, where K1=sum(Kj1) and K2=sum(Kj2), and the
%%%% columns are initially ordered based on a simple column stacking
%%%%
%%%% This function reorders them according to the J=(J1+1)*(J2+1) "levels", 
%%%% Where the order of the levels is done as follows:
%%%%
%%%%        Kj1(0)      Kj1(1)      Kj1(2)      ...         Kj1(J1)
%%%% Kj2(0)     1           2           5
%%%% Kj2(1)     3           4           6
%%%% Kj2(2)     7           8           9
%%%% ...        ...                             ...
%%%% Kj2(J2)
%%%%
%%%% That is, the first Kj1(0)*Kj2(0) coefficients are marked 1, then the
%%%% next Kj1(1)*Kj2(0) are 2, then Kj1(0)*Kj2(1) are 3, then Kj1(1)*Kj2(1)
%%%% are 4, Kj1(2)*Kj1(0) are 5, Kj1(2)*Kj1(1) are 6, Kj1(0)*Kj2(2) 7,
%%%% Kj1(1)*Kj2(2) 8, Kj1(2)*Kj2(2) 9, etc., going level by level, first
%%%% the column wavelets and then rows.   After repeating this for
%%%% min(J1,J2) levels, then the remaining just go down the appropriate
%%%% rows and columns to complete the transform for the remaining
%%%% levels
%%%% 
%%%% To use this, 
%%%%    1. Given vector d of wavelet coefficients in simple
%%%% vectorized order, let d=d(reorder) to reorder d according to levels as
%%%% described above.  This is the desired order for running FMM or doing 
%%%% regularization by wavelet levels
%%%%
%%%%    2. Given vector d of wavelet coefficients ordered according to
%%%%    wavelet levels as described above, let d(reorder)=d to reorder d
%%%%    in a vectorized order, which is what is needed to apply the
%%%%    rectangular 2d-idwt, which involves 1d-idwt to the columns and then
%%%%    to the rows.
%%%%

J1=length(Kj1);
J2=length(Kj2);
K1=sum(Kj1);
K2=sum(Kj2);
K1ctr_end=cumsum(Kj1);
K2ctr_end=cumsum(Kj2);
K1ctr_start=[1,K1ctr_end(1:(end-1))+1];
K2ctr_start=[1;K2ctr_end(1:(end-1))+1];
J=J1*J2;
Kj=zeros(1,J);
%%%% This may be a slow way to do this but will suffice for now 8/21/12
reorder=zeros(1,K1*K2);
    Dctr=1;
    jj=1;
    d=reshape(1:K1*K2,K2,K1);
    reorder(1:Kj1(1)*Kj2(1))=reshape(d(1:Kj2(1),1:Kj1(1)),1,Kj2(1)*Kj1(1));
    Kj(jj)=Kj1(1)*Kj2(1);
    Dctr=Dctr+Kj(jj);
    jj=jj+1;
    for j=2:min(J1,J2)
        for j2=1:(j-1)
            reorder(Dctr:(Dctr+Kj1(j)*Kj2(j2)-1))=reshape(d(K2ctr_start(j2):K2ctr_end(j2),K1ctr_start(j):K1ctr_end(j)),1,Kj1(j)*Kj2(j2));
            Kj(jj)=Kj1(j)*Kj2(j2);
            Dctr=Dctr+Kj(jj);
            jj=jj+1;
        end;
        for j1=1:j
            reorder(Dctr:(Dctr+Kj1(j1)*Kj2(j)-1))=reshape(d(K2ctr_start(j):K2ctr_end(j),K1ctr_start(j1):K1ctr_end(j1)),1,Kj1(j1)*Kj2(j));
            Kj(jj)=Kj1(j1)*Kj2(j);
            Dctr=Dctr+Kj(jj);
            jj=jj+1;
        end;
    end;
    if (J1>J2)
        for j1=(J2+1):J1
            for j2=1:J2
                reorder(Dctr:(Dctr+Kj1(j1)*Kj2(j2)-1))=reshape(d(K2ctr_start(j2):K2ctr_end(j2),K1ctr_start(j1):K1ctr_end(j1)),1,Kj1(j1)*Kj2(j2));
                Kj(jj)=Kj1(j1)*Kj2(j2);
                Dctr=Dctr+Kj(jj);
                jj=jj+1;
            end;
        end;
    elseif (J2>J1)
        for j2=(J1+1):J2
            for j1=1:J1 
                reorder(Dctr:(Dctr+Kj1(j1)*Kj2(j2)-1))=reshape(d(K2ctr_start(j2):K2ctr_end(j2),K1ctr_start(j1):K1ctr_end(j1)),1,Kj1(j1)*Kj2(j2));
                Kj(jj)=Kj1(j1)*Kj2(j2);
                Dctr=Dctr+Kj(jj);
                jj=jj+1;
            end;
        end;  
    end;



        

