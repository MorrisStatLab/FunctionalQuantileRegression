function [Vbetans_trans,W_trans,X1]=transform(q,xi,sigma,model,wavespecs,getW,getV)
% This function pre-shifts and scales D, pre-scales X to get the new
% model. 
T=wavespecs.T;

const1=(1-2*q)/q/(1-q);
const2=2/q/(1-q);

scale=const2*bsxfun(@times,xi,sigma');
Ng_scale=scale.^(-0.5);

X1=NaN(model.n,model.p,T);

if getW==1
   Y1=model.Y-const1*xi; 
   Y2=Ng_scale.*Y1;
   W_trans.XtX=NaN(model.p,model.p,T); 
   W_trans.XtD=NaN(model.p,T);   
   W_trans.DtD=Y2'*Y2;
else
   W_trans=[];
end

for t=1:T    
    Xnt=rowprod(Ng_scale(:,t),model.X);    
    X1(:,:,t)=Xnt;
    if getW==1
        W_trans.XtX(:,:,t)=Xnt'*Xnt;        
        W_trans.XtD(:,t)=Xnt'*Y2(:,t);
    end
end

if getV==1
    Vbetans_trans=NaN(model.p,T);
    for t=1:T
        Vbetans_trans(:,t)=diag(X1(:,:,t)'*X1(:,:,t)).^(-1);
    end        
else
    Vbetans_trans=[];
end
