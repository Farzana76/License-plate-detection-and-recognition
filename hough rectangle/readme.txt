This is not exactly the same code used in the paper "Rectangle Detection based 
on a Windowed Hough Transform", but works basically the same way. 


The main file is detect_rectangles_paa.m, and the usage is:

[rectangles, centers] = detect_rectangles_paa(img, RMin, RMax, amax);

where img is a binary (edge image), RMin and RMax are the internal and external
radii of the ring-like window, and amax is the largest expected value
for the minumim side of all rectangles (I guess this parameter is not
described in the original paper, but setting it to RMax should work).

Try the code with the .mat contained in the zipped file, using

[rectangles, centers] = detect_rectangles_paa(bw, RMin, RMax, RMax);

This should reproduce Figure 12 of the paper.