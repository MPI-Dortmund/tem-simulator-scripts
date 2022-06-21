// The purpose of this script is to create a 3D SPHIRE with image
// After using it you save it to disk and change the pixelsize to whatever you want using eman.
// Then you might can use it in TEMSimulator.
dimx=120;
dimy=120;
dimz=120;
value=1.0;
newImage("mask", "32-bit black", dimx, dimy, dimz);
setBatchMode(true);
cX = dimx/2;
cY = dimy/2;
cZ = dimz/2;
radius = 50

width = getWidth();
height = getHeight();


for (z = 1; z <= nSlices; z++) {
    setSlice(z);
    for (x=0;x < getWidth(); x++){
    	for (y = 0; y < getHeight(); y++) {
    		distance = sqrt(pow(cX-x,2)+pow(cY-y,2)+pow(cZ-z,2));
    		if (distance<radius) {
    			setPixel(x, y, value);
    		}
    	}
    }   
}
setBatchMode(false);