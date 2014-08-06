
# SPF-CellTracker

This is a software suite for tracking more than a hundred of cells from 4D live cell
imaging data with following characteristics:

1. gray scale image, 
2. globular-like cell shapes, 
3. strong covariation of cells to be tracked, 
4. no background image except background noise and salt- and-pepper noise. 

This software suite is composed of three executables: **convert**, **track**, and **view**.
The first and second software are implemented in C while the third is implemented in C++.
OpenCV 2.4+ is used for *convert*. OpenGL 3.3+ and freeglut are used for *view*.

## convert

The first software **convert** converts a set of 2D images that compose 4D live-cell
imaging data into a single file encoded as our original binary format. The software
*convert* utilizes OpenCV 2.4+ only for loading original image files, and thereby
all the image file formats to which OpenCV 2.4.9 supports can be converted. During
conversion, average subtraction and 3D median filter can be optionally applied for
each 3D image in order to remove background noise and salt and pepper noise. We
implemented a linear-time algorithm for computing median to make 3D median filter faster.

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

The parameters of *convert* are set in **conf-convert.txt** whose file format is as follows:

    prefix:image
    imtype:tif
    number:1,1,1
    imsize:512,256,20,1000
    sbmean:1
    median:1,1,1
    inform:1
    output:data.bin

- prefix: File name prefix of input 2D images. 
  - This parameter only accept any character sequence without white spaces.
- imtype: The image file format of input 2D image. 
  - Accessible file formats are those to which OpenCV 2.4+ supports. For example, *tif*, *png*, *jpg*, and so on.
- number: The type of file name and offset(s) for serial numbers. 
  - The first number specifies the type of input 2D file names.
  - The second and third number specify offset values for serial numbers.
  - Details are described later.
- imsize: The size of image. 
  - If `imsize:512,256,20,1000`, resolution of each 2D image, z-coordinate (depth), and the number of time points 
are 512x256, 20, and 1000, respectively.
- sbmean: Whether or not the average intensity of each 2D image is subtracted (1 or 0) from the image.
  - For removing the background noise whose level is dependent on the value of *z*-axis. 
- median: Specifies window size used for 3D median filter for removing salt and pepper noise. 
  - `median:l,m,n`--> The windows size become (2l+1, 2m+1,2n+1). 
  - `median:0,0,0`--> 3D median filter is skipped.
  - `median:1,1,0`--> 2D median filter with window size 3x3 for each 2D image.
- inform: Whether or not file names of input 2D images that are correctly read is printed. 
- output: Output file name.

The file names of input 2D images are customizable by the first four parameters *prefix*, *imtype*, *number*
and *imsize*.  Accessible file names of input 2D images are categorized into following two types:

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
In the second case, the first and third values must be 0. Serial numbers of *z* are incremented
in the order of *z*size\*t+z. These can start from any of non-negative integers specified by
the second value of `number`. 
 

## track
The second software **track** is the main software for automatic detection and tracking of
multiple cells based on the SPF algorithm from the converted 4D image file. For details of
the algorithm, please refer to the article titled *SPF-CellTracker: Multiple cell tracking with a
spatial particle filter* (in submission).

COMPILATION:

1. Get all the source codes registered in **tracker** directory of this
   repository and move them to the current directory.

2. Type `make` in the terminal window.

USAGE:

1. Modify the configuration file `conf-track.txt` according to your preference.

2. Type `./track` in the terminal window.

PARAMETERS:

The parameters of *track* can be set in **conf-track.txt**. The file format is as follows:

    seed:1
    cutoff:40
    dploop:10
    root:0
    N:1000
    alpha:0.7
    beta:0.05
    sgmt:0.6,0.6,0.03
    sgms:0.6,0.6,0.03
    wmax:1,1,1
    wlik:1,1,1
    objsize:7,7,2
    margin:3,3,0
    tune:4
    input:data.bin
    output:result.bin

- seed: The seed value for generating random numbers.
- cutoff: The cutoff value of fluorescent intensities (0 ... 255).
  - Required only for cell detection and counting.
- dploop: The number of iterations for DP-means algorithm.
  - DP-means algorithm is used for cell detection and counting. 
- root: The cell ID which is initially tracked for each time.
  - `root:0`--> The root cell ID is automatically defined. 
- N: The number of particles for tracking each cell.
- The balance parameter between covariation and relative positions among cells (0 ... 1.0).
- beta: Probability that a cell disappears for each time.
- sgmt: Standard deviation of system noise for the root cell.
  - The first, second, and third values corresponds to standard deviations according to *x*, *y*, and *z* coordinates, respectively.
- sgms: Standard deviation of system noise for non-root cells.
  - The first, second, and third values corresponds to standard deviations according to *x*, *y*, and *z* coordinates, respectively.
- wmax: Specifies window size required for cell detection.
  - Should be **narrower** than the size of cells.
  - `wmax:l,m,n`--> Window size becomes (2l+1, 2m+1, 2n+1).
- wlik: Specifies window size required for the computation of likelihoods.
  - `wlik:l,m,n`--> Window size becomes (2l+1, 2m+1, 2n+1).
- objsize: The size of cells.
  - `objsize:l,m,n`--> The cell size becomes (l, m, n).
  - The scale difference between *z* axis and *xy* plane is automatically controlled by this parameter. 
- margin: Defines area that avoids computation of likelihoods.
  - The parameter for tracking cells that are about to go out of the image space.
  - `margin:l,m,n`--> likelihoods are computed only if cells are within the region (l,m,n) voxels narrower than image space.
- tune: The parameter for the fine-tuning of tracked cell positions. 
  - Utilizes the result of the cell-detection by the *k*-means algorithm. 
  - Tracked cell positions are utilized as the initial parameters of the *k*-means algorithm.
  - `tune:4`--> If a tracked cell position is within 4 voxels from the position of the corresponding cell detected 
by the *k*-means algorithm, the cell position is modified to the position. 
- input: Input file name.
- output: Output file name.


## view

The third software **view** is to visualize the 4D image data with a
tracking result. This software allows us to zoom in, zoom out, and rotate
4D image during playing 4D image data. Currently, only a binary file for
Mac OS X is available for the software *view*. We will distribute source
codes in the near future.

INSTALL (Mac only):

1. Install the software **Xquartz**, an implementation of X window system 
   which runs on OS X. 

2. Get the binary file **view** in **bin/Mac** directory of this repository 
   and move it to the current directory.

USAGE:

1. Put a 4D image file converted by *convert* and a trajectory file output 
   by *tracker* in the current directory.

2. Type the following command:

  `./view input-file-name trajectory-file-name`





