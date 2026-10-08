//Be sure HiC image is selected!
//Should be 32 bit unsigned tiff. Scaled to 0.01 max.

//rotate and recolor
run("Rotate 90 Degrees Left");
run("Multiply...", "value=100000");
run("Enhance Contrast", "saturated=0.35");
run("Fire");
run("16-bit");

//copy the plot and name it
run("Duplicate...", "title=HiC_10kb.tif");

//find linear plot dimesnions
ImageSize=getHeight()
MinusImageSize=ImageSize*-1
PaddedImageSize=ImageSize+200

//make 100 pixel padding around plot to avoid edge effects
run("Canvas Size...", "width=PaddedImageSize height=PaddedImageSize position=Center");
run("Duplicate...", " ");
run("Translate...", "x=0 y=ImageSize interpolation=None");
selectImage("HiC_10kb.tif");
run("Duplicate...", " ");
run("Translate...", "x=0 y=MinusImageSize interpolation=None");
imageCalculator("Add", "HiC_10kb-1.tif","HiC_10kb-2.tif");
selectImage("HiC_10kb-2.tif");
close;
imageCalculator("Add", "HiC_10kb.tif","HiC_10kb-1.tif");
selectImage("HiC_10kb-1.tif");
close;
run("Duplicate...", " ");
run("Translate...", "x=ImageSize y=0 interpolation=None");
selectImage("HiC_10kb.tif");
run("Duplicate...", " ");
run("Translate...", "x=MinusImageSize y=0 interpolation=None");
imageCalculator("Add", "HiC_10kb-1.tif","HiC_10kb-2.tif");
selectImage("HiC_10kb-2.tif");
close;
imageCalculator("Add", "HiC_10kb.tif","HiC_10kb-1.tif");
selectImage("HiC_10kb-1.tif");
close;

//subtract background
run("Duplicate...", "title=HIC_10kb_background");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=10 create");
imageCalculator("Subtract create", "HiC_10kb.tif","HIC_10kb_background");
selectImage("Result of HiC_10kb.tif");
rename("HiC_10kb_RB.tif");

//blur image to avoid small peaks
run("Duplicate...", "title=HiC_10kb_RB_smoothed");
run("Median...", "radius=2");
run("Gaussian Blur...", "sigma=3");

//threshold to remove ridge
setAutoThreshold("Li dark no-reset");
run("Analyze Particles...", "size=2000-Infinity show=Masks display clear");
run("Gaussian Blur...", "sigma=1");
run("Max...", "value=1");
run("Multiply...", "value=256");
run("Invert");
run("Divide...", "value=256");
imageCalculator("Multiply create", "HiC_10kb_RB.tif","Mask of HiC_10kb_RB_smoothed");
selectImage("Result of HiC_10kb_RB.tif");
rename("HiC_10kb_RB_RidgeRemoved.tif");

//subtract background shot noise. Estimated as 2.5*mean
getStatistics(area, mean, min, max, std, histogram)
NewOffset=mean*2.5
run("Subtract...", "value=NewOffset");

//Create gaussian blurred image for printing at small resolution
run("Duplicate...", "title=HiC_blur_contrast.tif");
run("Gaussian Blur...", "sigma=1");
run("Canvas Size...", "width=ImageSize height=ImageSize position=Center");
setMinAndMax(0, 3000);
run("Invert LUT");

//Create profile of the peaks for comparison to ChIP and RNAseq. Save data from plot.
selectImage("HiC_10kb_RB_RidgeRemoved.tif");
run("Canvas Size...", "width=ImageSize height=ImageSize position=Center");
setMinAndMax(0, 3000);
run("Select All");
run("Plot Profile");


