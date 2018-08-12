deblur_adl code README
December 18, 2010

Code installation:
-------------------------
1. Download KSVDBox and OMPBox before the use of this code, the toolboxes
can be found on Ron Rubinstein's website 
http://www.cs.technion.ac.il/~ronrubin/software.html 
About how to install their code, please refer to their readme file.
2. Add the KSVDBox and OMPBox directory to the Matlab path.
3. Add the deblur_adl code installation directory to the Matlab path.

Quick start:
-------------------------
1. Enter "demo1" or "demo2" at the Matlab command window to run a demonstration

Notice and problems:
-------------------------
1. Make sure the patch number within the image window is great than the dictionary size 
(KSVD algorithm requires an over-complete dictionary)
2. The size of specified image window affects the computational speed, better not larger than
250x250.
3. The selection of specified image window affects the deblur result, better try several 
windows.
3. The convergence rate is quick, do not need to set up a large iteration number.
4. So far, the algorithm does not deal with large blur kernel well.