
# SPF-CellTracker

This is a software suite for tracking more than several hundreds of cells from 4D live cell
imaging data with following characteristics:

1. gray scale image, 
2. globular-like cell shapes, 
3. strong covariation of cells to be tracked, 
4. no background image except background noise and salt- and-pepper noise. 

<a href="http://www.youtube.com/watch?feature=player_embedded&v=-IFJG-LZmfs"
target="_blank"><img src="http://img.youtube.com/vi/-IFJG-LZmfs/0.jpg"
alt="Visualization of neurons in a nematode, C. elegans." width="300" height="240" border="10" /></a>
<a href="http://www.youtube.com/watch?feature=player_embedded&v=IdPZ3d_D-iM"
target="_blank"><img src="http://img.youtube.com/vi/IdPZ3d_D-iM/0.jpg"
alt="Visualization of neurons in a nematode, C. elegans." width="300" height="240" border="10" /></a>

Details of the algorithm is available [here](https://ieeexplore.ieee.org/document/8186251).
This software suite is composed of four executables: `spf-convert`, `spf-detect`, `spf-track`,
and `spf-view`. All of them are implemented in C/C++ with OpenMP, OpenCV, OpenGL, freeglut.

## Compilation:

1. Install XQuartz, GCC-mp-6, OpenCV, OpenGL, and freeglut.

2. Type `make` in the `src` directory. If you failed the compilation, try this:
  - Replace `gcc-mp-6` with `gcc` in `src/tracker/main.c`, which disables the use of OpenMP.

## Demo:

1. Compile source codes by following the above instruction.

2. Type `./spf-pipeline`.

## Pre-processing

The binary `spf-convert` converts 4D image sequence composed of multiple image files.
The 8-bit grayscale images can be converted. If your image sequence is composed of 16-bit grayscale
images, convert them into 8-bit images by using ImageMagick, for example.
The conversion is performed by the following Bash script with the ImageMagick command `convert`:

```
$TIME=200
$ZDEPTH=20
for t in `seq -f%04g 1 $TIME`; do 
  echo $t;
  for z in `seq -f%02g 1 $ZDEPTH`; do 
    convert yourimage_t${t}_z00${z}.tif -evaluate multiply 16 -depth 8 yourimage8bit_t${t}_z00${z}.tif    
  done;
done;
```

The parameters of `spf-convert` can be set in `conf-convert.txt` whose file format is as follows:

    prefix:image
    imtype:tif
    number:1,1,1
    imsize:512,256,20,1000
    sbmean:1
    median:1,1,1
    inform:1

- prefix: The prefix of file names for input 2D images. 
  - This parameter only accepts any character sequence without white spaces.
- imtype: The image file format of input 2D images. 
  - Any formats of 8-bit gray scale images are acceptable if they are supported by OpenCV.
    For example, *tif*, *png*, *jpg*, and so on.
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

The file names of input 2D images are can be customized by the first four parameters *prefix*, *imtype*, *number*
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

## Cell Detection

`spf-detect` is an executable for finding cells captured by a cell image with z-stacks.
It requires the result of `spf-convert`.
To start it up, follow the command instructed by `spf-pipeline.sh`.
Cells are detected by the following procedure:
1. Find local peaks of voxel values from the image with the Laplacian of Gaussian filter.
2. Cluster local peaks of voxel values within the maximum cell size.

The parameters of `spf-detect` can be set in ``setting/conf-detect.txt`` The file format is as follows:

    Intensity cutoff            10
    Max number of cells         512
    Peak bounding box           1,1,1
    Cell bounding box           10,10,1
    Lambda in dpmeans           5.0
    Number of loops in dpmeans  10
    Standard dev. in log        1.0

- Intensity cutoff: Voxel values which are less than this value are ignored.
- Max number of cells: The maximum number of cells to be reported.
- Peak bounding box: defines the window size of a local peak of voxel values.
  - `l,m,n` --> corresponds to the window size (2l+1, 2m+1, 2n+1).
  - Local peaks in the window are considered to be the same peak.
- Cell bounding box: defines the window size of the LoG filter, which can be interpreted as the maximum size of cells.
  - `l,m,n` --> corresponds to the window size (2l+1, 2m+1, 2n+1).
- Number of loops in dpmeans: The number of iterations for the DP-means algorithm, which
  estimates a cell position by clustering local peaks.
- Standard dev. in log: defines the blob size for the Laplacian of Gaussian filter.

## Cell tracking

`spf-track` is the executable that tracks cells captured by 4D image sequence.
It requires the results of `spf-convert` and `spf-detect`.
To start it up, follow the command instructed by `spf-pipeline.sh`.
The parameters of `spf-track` can be set in ``setting/conf-track.txt``. 
The file format is as follows:

    np:1000
    alpha:0.85
    sgmt:3.0,3.0,0.1
    sgms:1.5,1.0,0.1
    sgml:10
    wlik:5,5,1

- np: The number of particles.
- alpha: The balance parameter that controls importance of covariation and relative positions among cells (0 ... 1.0).
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

## Visualization

`spf-view` is an executable for visualizing 4D imaging sequence and tracking results.
It requires the result of `spf-convert`. The result of `spf-track` is optionally visualized.
To start it up, follow the command instructed by `spf-pipeline.sh`.

Usage:

-   `b`: Toggle display mode regarding trackers: all/one-by-one
-   `c`: Increase the cutoff of voxel intensities
-   `d`: Decrease the cutoff of voxel intensities
-   `e`: Enable/Disable the voxel aspect ratio along z-axis
-   `f`: Enable/Disable focus mode
-   `m`: Next tracker in focus mode
-   `n`: Previous tracker in focus mode
-   `o`: Draw/Delete raw intensity data
-   `q`: Exit viewer
-   `r`: Reset zoom
-   `s`: Start/Stop
-   `t`: Draw/Delete tracker positions
-   `w`: Reset rotation originated from keyboard operation
-   `x`: rotate around x-axis
-   `y`: rotate around y-axis
-   `z`: rotate around z-axis
-   `.`: Reset rotation originated from mouse operation
-   `,`: Reset translation originated from mouse operation
- `ESC`: Exit viewer

