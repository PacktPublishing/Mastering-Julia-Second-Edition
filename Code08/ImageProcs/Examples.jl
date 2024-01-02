### [ A tale of two bridges } ###

#=
  The image TB1.jpgis a color view of London's Tower Bridge.
  turns it into a greyscale image, blurs it with a convolution 
  using a 4x4 Gaussian kernel, and displays it, together with 
  the grayscale as a 1-row mosaic.
=#
using Images, ImageView, Colors 
using ImageFiltering, MosaicViews

DS = ENV["HOME"]*"MJ2/DataSources") 
img = load("$DS/Files/TB1.jpg")

# Display the image type
typeof(img) 

imgB = imfilter(img, Kernel.gaussian(4)) 
imgX = mosaicview(img, imgB, nrow=1, npad=20, fillvalue = Gray{N0f8}(1.0)) 
imshow(imgX);

# Edge detection of the Brooklyn bridge

using ImageEdgeDetection, ImageTransformations
img = load("NYB2.png") 

# Display the image type
typeof(img) 

imgR = resize(img, ratio=0.5); 
imgE = detect_edges(Gray.(imgR), Canny(spatial_scale=1.2))
imgX = mosaicview(imgR, imgE, nrow=1, npad=20, fillvalue = Gray{N0f8}(1.0)) 
imshow(imgX);
