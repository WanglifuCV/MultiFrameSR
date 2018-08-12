function b_im = Blur_M( im, B )

[h w ch]  =  size(im);

if ch==3
    m1  = im(:,:,1);
    m2  = im(:,:,2);
    m3  = im(:,:,3);    
    b_im(:,:,1)  =  reshape( B*m1(:), h, w );
    b_im(:,:,2)  =  reshape( B*m2(:), h, w );
    b_im(:,:,3)  =  reshape( B*m3(:), h, w );
else
    b_im  =  reshape( B*im(:), h, w );
end