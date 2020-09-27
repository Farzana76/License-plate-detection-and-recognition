I=imread('36dr.jpg');
I=rgb2gray(I);
I=imbinarize(I);
I=imcomplement(I);
imwrite(I,'1ch.bmp');