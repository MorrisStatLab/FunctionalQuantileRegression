function A=form_2d_dwt(theta,wlevels)

%%%% Form concatenated vector of wavelet coefficients for 2-d idwt.
%%%%    by entering in elements of theta as the diagonal elements of the
%%%%    2d-dwt object, representing the covariance structure of the 
%%%%    wavelet coefficients.
%%%%
%%%%    Inputs: 
%              cov = assumption on variance components
%                               (0=homoscesastic, 1=homoscedastic within levels, 2=heteroscedastic)
%%%%           theta: vector of length * (1, J, or K) containing
%%%%                    wavelet-space variance components
%%%%            wlevels: vector of number of coefficients per wavelet
%%%%                        level.
%%%%
%%%%    Outputs:    A = concatenated vector of wavelet coefficients, ready
%%%%                    for waverec2
%%%%

J=length(wlevels);
A=repmat(0,wlevels(1)^2+3*wlevels(2:end)'*wlevels(2:end),1);

kctr=1;
Arange=1:wlevels(1)^2;
A(Arange)=reshape(diag(theta(1:wlevels(1))),wlevels(1)^2,1);
kctr=kctr+wlevels(1);
Actr=Arange(end)+1;
for j=2:J
    Arange=(Actr+2*wlevels(j)^2):(Actr+3*wlevels(j)^2-1);
    A(Arange)=reshape(diag(theta(kctr:(kctr+wlevels(j)-1))),wlevels(j)^2,1);
    kctr=kctr+wlevels(j);
    Actr=Arange(end)+1;
end
