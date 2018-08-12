function [imout fft_h] = Generate_blur_image(im, Blur_type, Blur_par, nSig, seed)

[h w ch] =  size(im);

if Blur_type==1
    psf      =  ones(Blur_par)/(Blur_par*Blur_par);
    fft_h    =  zeros(h,w);
    t        =  floor( size(psf,1)/2 );
    fft_h(h/2+1-t:h/2+1+t,w/2+1-t:w/2+1+t)  = psf;
    fft_h    =  fft2(fftshift(fft_h));    
else
    % Generates the kernel in frequency domain
    
%     sigall = Blur_par^2*2;  
% 
%     %Gaussian Blurring Function
%     sigblurx=sigall;
%     sigblury=sigall;
% 
%     for i=1:h
%         for j=1:w
%              hhg(i,j)=exp(-(i-h/2-1).^2/sigblurx).*exp(-(j-w/2-1).^2/sigblury);
%         end
%     end
%     hh=hhg;
%     hh=hh/sum(sum(hh));
%     hf=fftshift(hh);
%     fft_h=fft2(hf);    

    psf=fspecial('gaussian', 25, Blur_par);
    fft_h    =  zeros(h,w);
    t        =  floor( size(psf,1)/2 );
    fft_h(h/2+1-t:h/2+1+t,w/2+1-t:w/2+1+t)  = psf;
    fft_h    =  fft2(fftshift(fft_h));        
    
%     u = meshgrid(-h/2:h/2-1,-w/2:w/2-1);
%     u = u/h;
%     v = u.';    
%     sigma_h = Blur_par;
%     H_f =  exp(-pi^2*sigma_h^2*(u.^2+v.^2));
%     H_f = fftshift(H_f);
%     psf   =  fftshift(real(ifft2(H_f)));
%     psf   =  psf(2:end,2:end);    
%     
%     fft_h    =  zeros(h,w);
%     t        =  floor( size(psf,1)/2 );
%     fft_h(h/2+1-t:h/2+1+t,w/2+1-t:w/2+1+t)  = psf;
%     fft_h    =  fft2(fftshift(fft_h));    
    
end

if ~exist('seed'),
    seed = 0;
end    
randn('seed',seed);


imout   =  zeros( size(im) );
if ch==3
    for i = 1 : 3
        im_f    =  fft2(im(:,:,i));
        z_f     =  fft_h.*im_f;
        z       =  real(ifft2(z_f));
        imout(:,:,i)    =  z;
    end
else
    im_f    =  fft2(im);
    z_f     =  fft_h.*im_f;
    z       =  real(ifft2(z_f));
    imout   =  z;
end
    
imout   =  imout + nSig*randn(size(im));
