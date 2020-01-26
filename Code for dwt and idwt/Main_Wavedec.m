% Example code to run the function DWT_rows

%First example DWT by rows
wavespecs.ndim=1;
%wavespecs.nlevels=5;          %%% applied to rows
wavespecs.boundary='per'; 
wavespecs.wavelet='db3';
load('train.mat');
X=train_Y;
wavespecs.compress=1;
[D1,wavespecs1]=DWT_rows(X,wavespecs);
Y1=IDWT_rows(D1,wavespecs1);
mean(mean(abs(Y1-X)))
[D1_compressed,wavespecs1]=wavelet_compress(D1,wavespecs1,0.999);
Y1C=IDWT_rows(D1_compressed,wavespecs1);
mean(mean(abs(Y1C-X)))

%Test for 2D rectangular wavelet decomposition
load('train.mat');
X=train_Y;
clear wavespecs
wavespecs.compress=1;
wavespecs.ndim=2;
wavespecs.rectangular=1; %%% rectangular 2d transform
%wavespecs.nlevels=[5,5];          %%% applied to [rows,cols]
wavespecs.boundary=['per';'sym']; 
wavespecs.wavelet=['db3';'db3'];
wavespecs.t=[120,120];
[D2r,wavespecs2r]=DWT_rows(X,wavespecs);
Yr=IDWT_rows(D2r,wavespecs2r);
mean(mean(abs(Yr-X)))
[D2r_compressed,wavespecs2r]=wavelet_compress(D2r,wavespecs2r,0.999);
YrC=IDWT_rows(D2r_compressed,wavespecs2r);
mean(mean(abs(YrC-X)))

%%% Now test 2d square transform
clear wavespecs
wavespecs.compress=1;
wavespecs.ndim=2;
wavespecs.rectangular=0; 
%wavespecs.nlevels=5;          
wavespecs.boundary='per'; 
wavespecs.wavelet='db3';
wavespecs.t=[120,120];
[D2s,wavespecs2s]=DWT_rows(X,wavespecs);
Ys=IDWT_rows(D2s,wavespecs2s);
mean(mean(abs(Ys-X)))
[D2s_compressed,wavespecs2s]=wavelet_compress(D2s,wavespecs2s,0.999);
YsC=IDWT_rows(D2s_compressed,wavespecs2s);
mean(mean(abs(YsC-X)))

%%% Now test 3d (square) transform
clear wavespecs
wavespecs.compress=1;
wavespecs.ndim=3;
wavespecs.rectangular=0; 
%wavespecs.nlevels=3;          
wavespecs.boundary='per'; 
wavespecs.wavelet='db2';
wavespecs.t=[20,24,30];
[D3,wavespecs3]=DWT_rows(X,wavespecs);
Y3=IDWT_rows(D3,wavespecs3);
mean(mean(abs(Y3-X)))
[D3_compressed,wavespecs3]=wavelet_compress(D3,wavespecs3,0.999);
Y3C=IDWT_rows(D3_compressed,wavespecs3);
mean(mean(abs(Y3C-X)))