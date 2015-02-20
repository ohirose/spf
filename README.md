
# SPF-CellTracker

This is a software suite for tracking more than a hundred of cells from 4D live cell
imaging data with following characteristics:

1. gray scale image, 
2. globular-like cell shapes, 
3. strong covariation of cells to be tracked, 
4. no background image except background noise and salt- and-pepper noise. 

<a href="http://www.youtube.com/watch?feature=player_embedded&v=JVTLNgwlkwg" 
target="_blank"><img src="http://img.youtube.com/vi/JVTLNgwlkwg/0.jpg" 
alt="Visualization of neurons in a nematode, C. elegans." width="300" height="240" border="10" /></a>
<a href="http://www.youtube.com/watch?feature=player_embedded&v=Bl71eWAzfDQ" 
target="_blank"><img src="http://img.youtube.com/vi/Bl71eWAzfDQ/0.jpg" 
alt="Visualization of neurons in a nematode, C. elegans." width="300" height="240" border="10" /></a>

This software suite is composed of three executables: **convert**, **track**, and **view**.
The first and second software are implemented in C while the third is implemented in C++.
OpenCV 2.4+ is used for *convert*. OpenGL 3.3+ and freeglut are used for *view*.

## convert

The first software **convert** converts a set of 2D images that compose 4D live-cell
imaging data into a single file encoded as our original binary format. The software
*convert* utilizes OpenCV 2.4+ only for loading original image files, and thereby
all the image file formats which OpenCV 2.4+ supports can be converted. 
During conversion, average subtraction for each 2D image and 3D median filter for
each set of 2D images that compose a 3D image can be optionally applied.
The average subtraction removes back ground noise whose level is strongly dependent on 
the depth of *z*-axis while the 3D median filter removes salt and pepper noise. 
We implemented a linear-time algorithm for computing median to make 3D median filter faster.

COMPILATION:

1. Install OpenCV 2.4+.

2. Get all the source codes registered in **converter** directory of this 
   repository and move them to the current directory.

3. Modify paths to OpenCV library and header files indicated by **LIBPATH** 
   and **INCLUDE** in **makefile** if needed. 

4. Type `make` in the terminal window.


USAGE:

1. Modify the configuration file **conf-convert.txt** according to your preference.

2. Type `./convert` in the terminal window.

PARAMETERS:

The parameters of *convert* can be set in **conf-convert.txt** whose file format is as follows:

    prefix:image
    imtype:tif
    number:1,1,1
    imsize:512,256,20,1000
    sbmean:1
    median:1,1,1
    inform:1
    output:data.bin

- prefix: The prefix of file names for input 2D images. 
  - This parameter only accepts any character sequence without white spaces.
- imtype: The image file format of input 2D images. 
  - Any formats that OpenCV 2.4+ supports are acceptable. For example, *tif*, *png*, *jpg*, and so on.
- number: Specifies a type of input file names and offsets for serial numbers. 
  - The first number specifies a type of input 2D file names.
  - The second and third number specify offset values for serial numbers.
  - Details are described later.
- imsize: The size of image. `imsize:512,256,20,1000`-->
  - Resolution of each 2D image is **256x512**. The first and second number specify width and height of each 2D image, respectively. 
  - The number of 2D images that compose a 3D image is 20.
  - The number of time points are 1000.
- sbmean: Whether or not the average intensity of each 2D image is subtracted from the image, aiming at
removing the background noise whose level is dependent on the value of *z*-axis. 
  - `sbmean:1`--> The average is subtracted.
  - `sbmean:0`--> The average is **NOT** subtracted.
- median: Specifies window size used for 3D median filter for removing salt and pepper noise. 
  - `median:l,m,n`--> The window size becomes (2l+1, 2m+1,2n+1). 
  - `median:0,0,0`--> 3D median filter is skipped.
  - `median:1,1,0`--> 2D median filter with window size 3x3 for each 2D image.
- inform: Whether or not file names of input 2D images correctly read is printed. 
  - `inform:1`--> File names are printed.
  - `inform:0`--> File names are **NOT** printed.
- output: Output file name.

The file names of input 2D images are customizable by the first four parameters *prefix*, *imtype*, *number*
and *imsize* of *conf-convert.txt*.  Accessible file names of input 2D images are categorized into following 
two types:

(i)

    image_t0001_z0001.tif
    image_t0001_z0002.tif
    ...
    image_t1000_z0019.tif
    image_t1000_z0020.tif

(ii)

    image1.tif
    image2.tif
    ...
    image19999.tif
    image20000.tif

The file name type can be switched by the first value of `number` (1 or 0).  In the first case, 
the first value of `number` must be 1. Symbols 't' and 'z' in the image file names represent 
time and value of *z*-axis, respectively. Serial numbers for time and *z*-value can start from any of non-negative
integers specified by the second and third values of `number`, respectively. The number of time points 
and the number of *z*-values are indicated by the third and fourth values of `imsize`. 
In the second case, the first and third values must be 0. Serial numbers are incremented
in the order of zsize\*t+z. These can start from any of non-negative integers specified by
the second value of `number`. 
 

## track
The second software **track** is the main software for automatic detection and tracking of
multiple cells based on the SPF algorithm from the converted 4D image file. For details of
the algorithm, please refer to the article titled *SPF-CellTracker: Multiple cell tracking with a
spatial particle filter* (in submission).

COMPILATION:

1. Get all the source codes registered in **tracker** directory of this
   repository and move them to your current directory.

2. If you are Linux user, replace **gcc-mp-4.8** with **gcc** in **makefile**.
   Type `make` in the terminal window.

USAGE:

1. Modify the configuration file `conf-track.txt` according to your preference.

2. Type `./track` in the terminal window.

PARAMETERS:

The parameters of *track* can be set in **conf-track.txt**. The file format is as follows:

    cutoff:18
    wmax:2,2,1
    lambda:6
    zscale:3
    alpha:0.85
    sgmt:3.0,3.0,0.1
    sgms:1.5,1.0,0.1
    sgml:10
    wlik:5,5,1
    input:data.bin
    output:result.bin

- cutoff: The cutoff value of fluorescent intensities (0 ... 255).
  - Required only for cell detection and counting.
- wmax: Specifies window size required for cell detection.
  - Should be **narrower** than the size of cells.
  - `wmax:l,m,n`--> Window size becomes (2l+1, 2m+1, 2n+1).
- lambda: The maximum radius of cells to be tracked.
- zscale: Ratio of the voxel step length along *z*-axis to that in *xy*-plane.
- alpha: The balance parameter that controlls importance of covariation and relative positions among cells (0 ... 1.0).
  - `alpha:1.0`--> Information of covariation among cells is only utilized.
  - `alpha:0.0`--> Information of initial relative positions among cells is only utilized.
  - `alpha:0.5`--> The both information is utilized with equal weight.
- sgmt: Standard deviation of system noise for the root cell.
  - The first, second, and third values corresponds to standard deviations according to *x*, *y*, and *z* coordinates, respectively.
- sgms: Standard deviation of system noise for non-root cells.
  - The first, second, and third values corresponds to standard deviations according to *x*, *y*, and *z* coordinates, respectively.
- sgml: Standard deviation of the likelihood function.
- wlik: Specifies window size required for the computation of likelihoods.
  - `wlik:l,m,n`--> Window size becomes (2l+1, 2m+1, 2n+1).
- input: Input file name.
- output: Output file name.


## view

The third software **view** is to visualize the 4D image data with a
tracking result. This software allows us to zoom in, zoom out, and rotate
4D image during playing 4D image data. Currently, only a binary file for
Mac OS X is available for the software *view*. We will distribute source
codes in the near future.

INSTALL (Mac only):

1. Install **XQuartz**, an implementation of X window system 
   which runs on OS X. 

2. Install **freeglut**, a library for OpenGL.

3. Get the binary file **view** in **bin/Mac** directory of this repository 
   and move it to the current directory.

USAGE:

1. Put a 4D image file converted by *convert* and a trajectory file output 
   by *tracker* in the current directory.

2. Run XQuartz and type the following command in your terminal window:

  `./view input-file-name trajectory-file-name`





