clear;
clc;

In_dir        =   'C:\TDDOWNLOAD\VOCtrainval_11-May-2009\VOCdevkit\VOC2009\JPEGImages';
Out_dir       =   'Data\Deblurring_Test_sets\Cropped_2';
fpath         =   fullfile(In_dir, '*.jpg');
im_dir        =   dir(fpath);
im_num        =   length(im_dir);

h0            =   234;
w0            =   280;

for i = 1:im_num
    
    I           =   double( imread(fullfile(In_dir, im_dir(i).name)) );
    [h w ch]    =   size(I);
    
    h1    =   h0;
    w1    =   w0;
    if  w<h
        h1   =   w0;
        w1   =   h0;
    end
    
    if  h<=h1 && w<=w1
        im   =  I;
    elseif h<h1
        sw   =  floor((w-w1)/2);
        im   =  I(:, sw:sw+w1-1, :); 
    elseif w<w1
        sh   =  floor((h-h1)/2);
        im   =  I(sh:sh+h1-1, :, :);  
    else
        sw   =  floor((w-w1)/2)+1;
        sh   =  floor((h-h1)/2)+1;
        im   =  I(sh:sh+h1-1, sw:sw+w1-1, :);  
    end            
    fname    =  strcat('p_', im_dir(i).name);
    
    imwrite(im./255, fullfile(Out_dir, fname));
end
