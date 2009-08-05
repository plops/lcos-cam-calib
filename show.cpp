/**
 * @file   1.cpp
 * @author Martin Kielhorn <kielhorn.martin@googlemail.com>
 * @date   Wed Jul 22 12:05:28 2009
 * 
 * @brief Display a grating on the phase only slm and capture the
 * diffraction pattern on camera.
 * was originally 0722/phase_scan
 * 
 */
#include <stdio.h>
#include <unistd.h>
#include <complex>
#include <iostream>
#include <typeinfo>

using namespace std;

const int do_bound_checks=1;


/** 
 * Macro for to access the img as a 1d array. 
 * Usage:
 * doall1d(img,img[i]=i)
 *
 * @param img Image that will be changed.
 * @param cmds 
 */
#define doall1d(img,cmds)			\
  for(int i=0;i<img.n;i++)			\
    do{cmds;}while(0)				

/** 
 * Macro for to access the img as a 2d array.
 * Usage:
 * doall(img,img(i,j)=i+j)
 * doall(img,double r=sqrt(i+0.3*j); img(i,j)=r)
 * 
 * @param img Image that will be changed.
 * @param cmds 
 */
#define doall(img,cmds)				\
  for(int i=0;i<img.w;i++)			\
    for(int j=0;j<img.h;j++)			\
      do{cmds;}while(0)				


/** 
 * Macro to access 2d image. It computes floating point coordinates x
 * and y with a range [-0.5,0.5).
 * Usage:
 * doalld(img,r=hypot(x,y); img(i,j)=j0(r))
 * 
 * @param img 
 * @param cmds 
 */
#define doalld(img,cmds)			\
  do{						\
  double x=-0.5,y,dx=1./img.w,dy=1./img.h;	\
  int j;					\
  for(int i=0;i<img.w;x+=dx,i++)		\
    for(j=0,y=-0.5;j<img.h;y+=dy,j++){		\
      cmds;					\
    }						\
  }while(0)				


class NoMemory{};
class NotInImage{};
class NotEvenDimension{};

template<class T> class Image {
public:
  const int w,h,n;
  T*data;
  unsigned int id;		/**< Space to store arbitrary
				   information. It is initialized to
				   0. Later I use this to store the
				   OpenGL texture object. */
  
  /** 
   * Construct a new WxH image. Space for data is allocated.
   * 
   * @param W 
   * @param H 
   */
  Image(const int W,const int H):w(W),h(H),n(W*H),data(0),id(0){
    data=new T[w*h];
  }

  /** 
   * Copy constructor. Its actually not copying data.
   * 
   * @param img 
   */
  Image(const Image<T>&img){
    this.data=0;
    *this=img;
  }

  ~Image(){
    delete [] data;
  }
  
  /** 
   * Access 2d image at column i and row j.
   * 
   * @param i column
   * @param j row
   * 
   * @return 
   */
  inline T&operator()(int i,int j){
    if(do_bound_checks)
      in_image_check(i,j);
    return data[i+w*j];
  }
  /** 
   * Consecutive (1D) access of image.
   * 
   * @param i 
   * 
   * @return 
   */
  T&operator[](int i) throw(NotInImage){
    if(do_bound_checks)
      if(i<0 || i>w*h)
        throw NotInImage();
    return data[i];
  }

  /** 
   * Assignmentoperator. Copy the data from input image.
   * 
   * @param img Input image
   */
  Image<T>&
  operator=(const Image<T>&img){
    if(this==&img)  // in case pic1=pic1 we don't have to copy anything
      return;
    if(data)
      delete [] data;
    w=img.w;
    h=img.h;
    n=w*h;
    id=0;
    data=new T[n];
    if(!data)
      throw NoMemory();
    else
      memcpy(data,img.data,n);
    return *this; // so that a=b=c works
  }
  // I don't know how to declare this function correctly Outside it is
  // templated but here it has to be T.  
  // friend ostream& operator<<(ostream&os,const Image<T>&img);

  /** 
   * Hackish way to get a printable representation of the template type.
   * like "float" or "double".
   * 
   * @return 
   */
  string print_T() {
    string str=__PRETTY_FUNCTION__; 
    // str looks like that: std::string Image<T>::print_T() [with T = double]
    size_t 
      pos=str.find(" = "),
      start=pos+3,
      end=str.find(']',start),
      diff=end-start;
    return str.substr(start,diff);
  }

private:
  void in_image_check(int i,int j) throw (NotInImage) {
    if(i<0 || i>=w || j<0 || j>=h)
      throw NotInImage();
  }
};

/** 
 * Add two images pixelwise. Creates a new image.  Use the bigger type
 * in the front as it will be the result type.
 * 
 * @param a 
 * @param b 
 * 
 * @return 
 */
template<class T,class S>Image<T>
operator+(const Image<T>&a,const Image<S>&b)
{
  if(a.w!=b.w || a.h!=b.h)
    throw NotInImage();
  Image<T> result(a.w,a.h);
  doall1d(result,result[i]=a[i]+b[i]);
  return result;
}

/** 
 * Subtract two images pixelwise. Creates a new image.  Use the bigger
 * type in the front. It will be the type of the result.
 * 
 * @param a 
 * @param b 
 * 
 * @return 
 */
template<class T,class S>Image<T>
operator-(const Image<T>&a,const Image<S>&b)
{
  if(a.w!=b.w || a.h!=b.h)
    throw NotInImage();
  Image<T> result(a.w,a.h);
  doall1d(result,result[i]=a[i]-b[i]);
  return result;
}

/** 
 * Function so that a textual representation of the image can be
 * printed out with cout.
 * 
 * @param os 
 * @param img  
 * 
 * @return Changed stream
 */
template<class T> ostream& 
operator<<(ostream&os,Image<T>&img){
  os<<"(Image "<<img.w<<" "<<img.h
    <<" ("<<img.print_T()<<") "<<img.id<<")";
  return os;
}


/** 
 * Find extremal pixel values in Image a.
 * 
 * @param a Image
 * @return extremal values in min and max 
 */
template<class T>void
extrema(Image<T>&a,T&min,T&max)
{
  min=max=a[0];
  for(int i=1;i<a.w*a.h;i++){
    T v=a[i];
    if(v>max)
      max=v;
    if(v<min)
      min=v;
  }
}

/** 
 * Run the function f for each pixel in the image src and write output
 * into dst.
 * 
 * @param src input image
 * @return dst output image
 * @param f function taking a value from the src image and returning a
 * value of for the dst image.
 */
template<class srcT,class dstT>void
transfer(Image<srcT>&src,
         Image<dstT>&dst,
         dstT(*f)(const srcT&))
{
  for(int i=0;i<src.n;i++)
    dst[i]=f(src[i]);
}

/** 
 * Reduce the function f over each of the columns of the image
 * src. Store the results in the (src.w)x1 image dst. The dst type
 * might be bigger than the source type.
 * 
 * @param src Input image
 * @param dst Row of values with width src.w.
 * @param f Function taking two values and generating a new one.
 */
template<class srcT,class dstT>void
reduce_columns(Image<srcT>&src,
	       Image<dstT>&dst,
	       dstT(*f)(dstT,dstT))
{
  if(src.w>dst.w)
    throw NotInImage();
  for(int i=0;i<src.w;i++){
    dstT sum=src(i,0);
    for(int j=1;j<src.h;j++)
      sum=f(sum,src(i,j));
    dst[i]=sum;
  }
}

// src.w,src.h -> project with f over src.h -> dst.w

template<class srcT,class dstT>void
reduce_rows(Image<srcT>&src,
	    Image<dstT>&dst,
	    dstT(*f)(dstT,dstT))
{
  if(src.h>dst.w)
    throw NotInImage();
  for(int i=0;i<src.h;i++){
    dstT sum=src(0,i);
    for(int j=1;j<src.w;j++)
      sum=f(sum,src(j,i));
    dst[i]=sum;
  }
}

/** 
 * Output unsigned char image one to one into pgm file.
 * 
 * @param a Image
 * @param filename 
 */
void
write_pgm(Image<unsigned char>&a,char*filename)
{
  FILE*f=fopen(filename,"w");
  fprintf(f,"P5\n%d %d\n255\n",a.w,a.h);
  if(a.h!=(int)fwrite(a.data,a.w,a.h,f))
    cerr<<"error writing unsigned char " 
	<< a.w << "x" << a.h << " image to "
	<<filename<<endl;
  fclose(f);
}

/** 
 * Write out a (real) image into filename scaling by min and max
 * 
 * @param a Image
 * @param filename 
 * @param min 
 * @param max 
 */
template<class T> void
write_pgm(Image<T>&a,char*filename,T min, T max)
{
  Image<unsigned char>buf(a.w,a.h);
  
  double s=255./(max-min); // if I use T instead of double, unsigned
			   // short doesn't work
  doall1d(buf,buf[i]=(a[i]-min)*s);

  write_pgm(buf,filename);
}

/** 
 * Write out a (real) image into filename so that the its minimal
 * pixel value will be 0 and its maximal pixel value 255.
 * 
 * @param a Image
 * @param filename 
 */
template<class T> void
write_pgm(Image<T>&a,char*filename)
{
  T min,max;
  extrema(a,min,max);
  write_pgm(a,filename,min,max);
}

/** 
 * Output complex image as absolute magnitude. Minimum value goes to 0
 * and maximum pixel value goes to 255.
 * 
 * @param a 
 * @param filename 
 */
void
write_pgm(Image<complex<double> >&a,char*filename)
{
  // temporary array with magnitude
  Image<double> scal(a.w,a.h);
  transfer(a,scal,abs);
  write_pgm(scal,filename);
}


#include <SDL/SDL.h>
#include <SDL/SDL_opengl.h>

/** 
 * Setup OpenGL coordinates to [0,0] in lower left and [w,h] in upper
 * right (setting for 2d drawing).
 * 
 * @param w width
 * @param h height
 */
void
reshape(int w,int h)
{
  glViewport(0,0,w,h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0,w,0,h,-1,1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

class Graphics{
public:
  int w,h;
  SDL_Surface*screen;
  /** 
   * Open an SDL OpenGL window.
   * 
   * @param W width
   * @param H height
   */
  Graphics(int W,int H):w(W),h(H){
    //SDL_putenv("SDL_VIDEO_X11_WMCLASS=littlegptracker"); // for eeepc
    SDL_putenv("SDL_VIDEO_WINDOW_POS=668,-31"); //1280-512-100
    SDL_Init(SDL_INIT_VIDEO);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER,1);
    screen=
      SDL_SetVideoMode(w,h,24,SDL_OPENGL);
    reshape(w,h);
  };
  ~Graphics(){
    SDL_FreeSurface(screen);
    SDL_Quit();
  };
};



const int target=GL_TEXTURE_RECTANGLE_NV; /**< I only use this as a
					     texture target */

void
prepare_load(unsigned int&id)
{
  if(id)
    glDeleteTextures(1,&id);
 
  glGenTextures(1,&id);
  glBindTexture(target,id);
  glEnable(target);
  //glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  int mode=GL_REPEAT; // GL_CLAMP_TO_EDGE
  glTexParameteri(target,GL_TEXTURE_WRAP_S,mode);
  glTexParameteri(target,GL_TEXTURE_WRAP_T,mode);
  mode=GL_LINEAR;
  glTexParameteri( target,
                   GL_TEXTURE_MIN_FILTER, 
                   mode);
  glTexParameteri( target,
                   GL_TEXTURE_MAG_FILTER,
                   mode );
  //  glPixelStorei(GL_UNPACK_ALIGNMENT,1);
}  

/** 
 * Upload unsigned char image as an OpenGL texture.
 * 
 * @param img 
 */
void
load(Image<unsigned char>&img)
{
  prepare_load(img.id);
  glTexImage2D(target,
               0, GL_LUMINANCE,
	       img.w,img.h,
               0,
               GL_LUMINANCE,
               GL_UNSIGNED_BYTE,
               img.data);
}

/** 
 * Destroy the OpenGL texture in the image and set the id to zero
 * again. This function should be used to mark an image as changed
 * after it has been modified.
 * 
 * @param img 
 */
template<class T>void
unload(Image<T>&img)
{
  glDeleteTextures(1,&img.id);
  img.id=0;
}

/** 
 * Upload unsigned short image as an OpenGL texture.
 * 
 * @param img 
 */
void
load(Image<unsigned short>&img)
{
  prepare_load(img.id);
  glTexImage2D(target,
               0, GL_LUMINANCE,
	       img.w,img.h,
               0,
               GL_LUMINANCE,
               GL_UNSIGNED_SHORT,
               img.data);
}

/** 
 * Draw image as an OpenGL texture. The origin of the texture
 * coordinates is in the upper left. If the given dimensions [w,h] are
 * bigger than the dimensions [img.w,img.h] the texture pattern will
 * be repeated to fill a bigger quad. E.g. the following code will
 * expand 320x1 image in y direction to display a 240 pixel high
 * grating:
 *
 * Image<unsigned char> img(320,1);
 * for(int i=0;i<img.w;i++)
 *   img(i,0)=(i%7)*(255/7);
 * draw(img,img.w,240);
 * 
 * @param img 
 * @param w 
 * @param h 
 */
template<class T>void
draw(Image<T>&img, int w, int h)
{
  if(img.id)
    glBindTexture(target,img.id);
  else
    load(img);

  int wt=w,ht=h;

  glEnable(target);
  glColor3d(1,1,1);
  
  glBegin(GL_QUADS);
  glTexCoord2i(0,ht);
  glVertex2i(0,0);

  glTexCoord2i(wt,ht);
  glVertex2i(w,0);

  glTexCoord2i(wt,0);
  glVertex2i(w,h);

  glTexCoord2i(0,0);
  glVertex2i(0,h);
  glEnd();
  glDisable(target);
}

/** 
 * Draw Image as an OpenGL texture.
 * 
 * @param img 
 */
template<class T>void
draw(Image<T>&img)
{
  draw(img,img.w,img.h);
}

#include <pvcam/master.h>
#include <pvcam/pvcam.h>

class PVCAM_Error{};

class Camera {
public:
  int16 hCam;
  Image<unsigned short>*frame;
  /** 
   * Open the Cascade II camera. It will deliver 512x512 images.
   * 
   */
  Camera(){
    frame=new Image<unsigned short>(512,512);

    // read ./src/hardware/cascadeII/pvcam26_57-062-00.pdf
    pl_pvcam_init();
    char cam_name[CAM_NAME_LEN];
    pl_cam_get_name(0,cam_name);
    pl_cam_open(cam_name,&hCam,OPEN_EXCLUSIVE);
  
    unsigned int a=
      //READOUT_PORT_MULT_GAIN;
      READOUT_PORT_NORMAL;
    pl_set_param(hCam,PARAM_READOUT_PORT,(void*)&a);
    //a=3000;
    //pl_set_param(hCam,PARAM_GAIN_MULT_FACTOR,(void*)&a);
    a=2;
    pl_set_param(hCam,PARAM_SPDTAB_INDEX,(void*)&a);
    a=CLEAR_PRE_SEQUENCE;
    pl_set_param(hCam,PARAM_CLEAR_MODE,(void*)&a);
    //uns16 exptime;
    //pl_get_param(hCam,ATTR_CURRENT,PARAM_EXP_TIME,(void*)&exptime);
    //printf("exposure time=%d\n",exptime);
    //pl_set_param(hCam,PARAM_EXP_TIME,(void*)&exptime);
    pl_exp_init_seq();
    const rgn_type region={0,511,1,0,511,1};
    uns32 size;
    const uns32 exptime=20; // depending on time resolution // perhaps millisecond
    // PARAM_EXP_RES and PARAM_EXP_RES_INDEX
    const int triggermode=TIMED_MODE;//TRIGGER_FIRST_MODE;
    pl_exp_setup_seq(hCam,1,1,&region,triggermode,exptime,&size);
  }
  ~Camera(){
    pl_exp_uninit_seq();
    pl_cam_close(hCam);
    pl_pvcam_uninit();
    delete frame;
  }
  /** 
   * Capture one image into frame.
   *
   * @return 512x512 Image that contains a pointer to the data.
   */
  Image<uns16>&
  capture(){
    pl_exp_start_seq(hCam,frame->data);
    int16 status;
    uns32 not_needed;
    while(pl_exp_check_status(hCam,&status,&not_needed) &&
	  (status!=READOUT_COMPLETE && status!=READOUT_FAILED))
      usleep(10000);
    if(status==READOUT_FAILED){
      cerr<<"data collection error: "<<pl_error_code()<<endl;
      throw PVCAM_Error();
    }
    pl_exp_finish_seq(hCam,frame->data,0); // decode data (would be
				     // important for multiple
				     // regions)
    return *frame;
  }
  /** 
   * Read out the current temperature of the camera sensor.
   * 
   * 
   * @return Temperature in degree celsius.
   */
  double
  temperature()
  {
    char strErr[ERROR_MSG_LEN];
    int16 temp = 0;
    if(pl_get_param(hCam,PARAM_TEMP,ATTR_CURRENT,&temp)==PV_FAIL) {
      pl_error_message(pl_error_code(),strErr);
      printf("error: %s\n",strErr);
    }
    return temp/100.;
  }
};

#define HAVE_INT32 // otherwise warning in master.h
#include <tiff.h>
#include <tiffio.h>

/** 
 * Output unsigned short 16 image into filename as TIFF. 
 * 
 * @param img 
 * @param filename 
 */
void
write_tiff(Image<uns16>&img,char*filename)
{
  TIFF*tif=TIFFOpen(filename,"w");
  
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, (uint32) img.w);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, (uint32) img.h); // height of the image

  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
  
  TIFFSetField(tif, TIFFTAG_COMPRESSION, 1);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, 1);
  TIFFSetField(tif, TIFFTAG_IMAGEDESCRIPTION, "phase_scan");
  TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
  int sampleperpixel=1;
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, sampleperpixel);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, 1);
  TIFFSetField(tif, TIFFTAG_DATETIME, "22/07/2009 12:00:00");
  TIFFSetField(tif, TIFFTAG_ARTIST, "Martin Kielhorn");

  // set strip size to be the size of one row of pixels
  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP,
	       TIFFDefaultStripSize(tif, img.w*sampleperpixel));

  // use most basic data storing method of the library, store one line
  // (row) at a time

  for(int row=0;row<img.h;row++)
    TIFFWriteScanline(tif,(tdata_t)&img(0,row),row,0);

  TIFFClose(tif);
}

unsigned short
max(unsigned short a,unsigned short b)
{
  if(a>b)
    return a;
  else
    return b;
}

/** 
 * Use region to store region of interests or with copy to crop
 * images.
 * 
 * @param cx Center of the rectangular region
 * @param cy 
 * @param W Width
 * @param H Height
 */
class Region{
public:
  int w,h,x,y;
  Region(double cx,double cy,int W,int H)
    :w(W),h(H),x(int(cx-W/2.)),y(int(cy-H/2.)){}
};

/** 
 * Function so that a textual representation of the region can be
 * printed out with cout.
 * 
 * @param os 
 * @param region
 * 
 * @return Changed stream
 */
ostream& 
operator<<(ostream&os,Region r){
  os<<"(Region "<<r.x<<" "<<r.y<<" "<<r.w<<" "<<r.h<<")";
  return os;
}


/** 
 * Draw the outline of the region.
 * 
 * @param r region
 */
void
draw(Region r)
{
  glPushMatrix();
  glTranslated(0,512,0);
 
  glScaled(1,-1,0); // needs to be upside down and shifted by 512
		    // because the [0,0] is in upper left of image
  glBegin(GL_LINE_LOOP);
  glVertex2i(r.x,r.y);
  glVertex2i(r.x+r.w,r.y);
  glVertex2i(r.x+r.w,r.y+r.h);
  glVertex2i(r.x,r.y+r.h);
  glEnd();
  glPopMatrix();
}

/** 
 * fill the 1d texture with a rect grating of given periodicity with
 * gray values either beeing min_val or max_val.
 * 
 * @param img 
 * @param period 
 * @param min_val 
 * @param max_val 
 */
void
fill_grating(Image<unsigned char>&img,
	     int period,int min_val,int max_val)
{
  for(int i=0;i<img.w;i++)
    img[i]=((max_val-min_val)*(i%period>period/2))+min_val;
  unload(img);
}

/** 
 * Copy a region from the input image into a another image. The
 * destination image has to have the right size.
 * 
 * @param src Input image
 * @param dst Input image
 * @param r Region to copy out of the src image
 * 
 * @return copied part of the src image
 */
template<class T>void
copy(Image<T>&src,Image<T>&dst,Region&r)
{
  for(int i=r.x;i<r.x+r.w;i++)
    for(int j=r.y;j<r.y+r.h;j++)
      dst(i-r.x,j-r.y)=src(i,j);
}

#define keydown(key,cmds)						\
      do{								\
	SDL_Event e;							\
	while(SDL_PollEvent(&e)){					\
	  if(e.type==SDL_KEYDOWN && e.key.keysym.sym==key){		\
	    cmds;							\
	  }								\
	}								\
      }while(0)								\

#include <vector>

vector<pair<int,int> >
find_max(Image<unsigned short>&img){
  if(img.h!=1) // this function is only defined for 1d vectors
    throw NotInImage();
  
  vector<pair<int,int> > peaks;
  for(int i=7;i<img.n-7;i++){
    unsigned short v=img[i];
    if(//(img[i-7] < v) &&
       (img[i-3] < v) &&
       (img[i-1] < v)
       && (img[i+1] < v) 
       && (img[i+3] < v)
       // && (img[i+7] < v)
)
      peaks.push_back(pair<int,int>(i,v));
  }
  return peaks;
}

double
integrate(Image<unsigned short>&img,
	  Image<double>&dark,
	  Region&r)
{
  double sum=0;
  for(int i=r.x;i<r.x+r.w;i++)
    for(int j=r.y;j<r.y+r.w;j++)
      sum+=img(i,j)-dark(i,j);
  return sum;
}

class Quad {
public:
  int a[4];
  Quad(int*q){
    for(int i=0;i<4;i++)
      a[i]=q[i];
  }
};

// hough transform stuff
enum{
    NTHETA=136,
    NRHO=138
};

unsigned char hough_buf[NTHETA*NRHO];
unsigned int hough_hist[NTHETA*NRHO]; 

const double 
theta_min=-M_PI/2,
  theta_max=M_PI/2,
  rho_max=M_SQRT2*512,
  rho_min=-rho_max;

// linear interpolation
inline double
stretch(int i,int n,double min,double max)
{
  double t=i*1./n;
  return (1-t)*min+t*max;
}

inline double
stretch_inclusive(int i,int n,double min,double max)
{
  double t=i*1./(n-1);
  return (1-t)*min+t*max;
}

inline double
theta(int i)
{
  return stretch(i,NTHETA,theta_min,theta_max);
}

inline double
rho(int i)
{
  return stretch_inclusive(i,NRHO,rho_min,rho_max);
}

#include <string.h>
double cos_tab[NTHETA],sin_tab[NTHETA];
void clear_hough()
{
  memset(hough_hist,0,sizeof(hough_hist));
}
void init_hough()
{
  clear_hough();
  for(int i=0;i<NTHETA;i++){
    double t=theta(i);
    cos_tab[i]=cos(t);
    sin_tab[i]=sin(t);
  }
}

void
insert_hough(int x,int y,unsigned int*hist)
{
  const static double sfrho=NRHO*1./(rho_max-rho_min);
  for(int i=0;i<NTHETA;i++){
    double rho=x*cos_tab[i]+y*sin_tab[i];
    hist[i+NTHETA*((int)((rho-rho_min)*sfrho))]++;
  }
}

// least squares line fitting

// collect points that are in the bin [maxtheta,maxrho] of the hough transform
void
extract_line_points(int maxtheta,int maxrho,vector<int>&x,vector<int>&y,
		    Image<unsigned char>&p)
{
  const static double sfrho=NRHO*1./(rho_max-rho_min);
  for(int j=0;j<p.h;j++)
    for(int i=0;i<p.w;i++){
      if(p(i,j)==0){
        double rho=j*cos_tab[maxtheta]+i*sin_tab[maxtheta];
        int irho=(int)((rho-rho_min)*sfrho);
        if(maxrho==irho){
          x.push_back(i);
          y.push_back(j);
        }
      }
    }
}
double
calc_avg(vector<int>x)
{
  double sum=0;
  for(unsigned int i=0;i<x.size();i++)
    sum+=x[i];
  return sum/x.size();
}

double
calc_x2(vector<int>x)
{
  double sum=0;
  for(unsigned int i=0;i<x.size();i++)
    sum+=x[i]*x[i];
  return sum;
}

double
calc_xy(vector<int>x,vector<int>y)
{
  double sum=0;
  for(unsigned int i=0;i<x.size();i++)
    sum+=x[i]*y[i];
  return sum;
}

class Error{};
// y=a+b*x  -> find a and b
double
fit_line(vector<int>x,vector<int>y,double&a,double&b)
{
  if(x.size()!=y.size() || x.size()==0)
    throw Error();
  int n=x.size();
  double 
    xbar=calc_avg(x),x2=calc_x2(x),
    ybar=calc_avg(y),y2=calc_x2(y),
    xy=calc_xy(x,y);
  double
    ssxx=x2-n*xbar*xbar,
    ssyy=y2-n*ybar*ybar,
    ssxy=xy-n*xbar*ybar;
  /*  cout << "xbar=" << xbar << " ybar="
       << ybar << " x2="
       << x2 << " y2="
       << y2 << " xy="
       << xy << " ssxx="
       << ssxx << " ssxy="
       << ssxy << " "<<endl;
  */
  b=ssxy/ssxx;
  //cout<<"b="<<b<<endl;
  a=ybar-b*xbar;
  //cout<<"a="<<a<<endl;
  double r2=ssxy*ssxy/(ssxx*ssyy);
  /* cout<<"correlation coeff r^2="
      <<r2
      <<" proportion of ssyy which is accounted for by the refgression"
      <<endl;
  */
return r2;
}


// x=a+b*y
void
draw_line_vert(double a,double b,Image<unsigned char>&p)
{
  for(int i=0;i<p.h;i++){ 
    double x=a+b*i;
     int ix=(int)x;
    if(ix>=3 && ix<p.w-3){
      p(ix-3,i)=128;
      p(ix+3,i)=128;
    } 
  }  
}

// y=a+b*x
void
draw_line_hori(double a,double b,Image<unsigned char>&p)
{
  for(int i=0;i<p.w;i++){ 
    double y=a+b*i;
    int iy=(int)y;
    if(iy>=3 && iy<p.h-3){
      p(i,iy-3)=128; 
      p(i,iy+3)=128; 
    }
  }  
}

// y=a+b*x
void draw_hor(pair<double,double>p)
{
  double a=p.first,b=p.second,x=0;
  glVertex2d(x,512-(a+b*x)); // in y 512-.. because of coordinate systems!
  x=512;
  glVertex2d(x,512-(a+b*x));
}
// x=a+b*y
void draw_ver(pair<double,double>p)
{
  double a=p.first,b=p.second,y=0;
  glVertex2d(a+b*y,512-y); // in y 512-.. because of coordinate systems!
  y=512;
  glVertex2d(a+b*y,512-y);
}

struct CoordMap {
  double x1,y1,x2,y2;
};

int
main()
{
  const int n=512;
  Graphics g(n+100+1920,1080); // open window
  Camera cam;

  vector<CoordMap>coord;
  FILE*f=fopen("scan20_thin/intersections.dat","r");
  // fill coord with contents from file
  double a,b,c,d;
  while(EOF!=fscanf(f,"%lg %lg %lg %lg",&a,&b,&c,&d))
    coord.push_back((CoordMap){b,a,c,d});
  fclose(f);

  Image<unsigned char> buf(n,n);

  Image<unsigned char> grating(100+1920,1);
  fill_grating(grating,20,0,255);
  
  for(int cnt=0;;cnt++){    
    cam.capture();
  
    glClear(GL_COLOR_BUFFER_BIT);
    glPushMatrix();
    
    // draw camera image in top left
    glPushMatrix();
    
    glTranslated(0,g.h-cam.frame->h,0);
    
    //unload(*cam.frame); draw(*cam.frame); //force update of texture
    double max_val=(*cam.frame)(0,255);
    for(int i=1;i<n;i++){
      double v=(*cam.frame)(i,255);
      if(v>max_val)
	max_val=v;
    }
    max_val*=1.6;
    unsigned int dark=1700;
    double sfbuf=255./(max_val-dark);
    doall1d(buf,buf[i]=(unsigned char)(((*cam.frame)[i]-dark)*sfbuf));
    unload(buf);
    draw(buf);
    {
    int cw=15,l=0;
    for(int j=0;l<coord.size();j++){
      glBegin(GL_LINES);
      for(int i=0;i<cw;i++){
	glVertex2d(coord[i+cw*j].x2,
		   512-coord[i+cw*j].y2);
	l++;
      }
      glEnd();
    }
    }
    // draw horizontal cross section from center line
    double graph_h=200;
    glTranslated(0,-graph_h,0);
    glBegin(GL_LINE_STRIP);
    for(int i=0;i<n;i++)
      glVertex2d(i,(*cam.frame)(i,255)*graph_h/(1<<16));
    glEnd();
    glColor3d(1,.2,.2);
    glBegin(GL_LINES);
    const int nstrips=10;
    for(int i=0;i<nstrips;i++){
      double h=i*graph_h/(nstrips-1);
      glVertex2d(0,h);glVertex2d(512,h);
    }
    glColor3d(1,.8,0);
    // this is the maximum which is displayed right now
    double y=max_val*graph_h/(1<<16);
    glVertex2d(0,y);
    glVertex2d(512,y);
    y=dark*graph_h/(1<<16);
    glVertex2d(0,y);
    glVertex2d(512,y);
    glEnd();
    glColor3d(1,1,1);
    
    glPopMatrix();
    
    glTranslated(512,0,0);
    //draw(grating,grating.w,g.h);
    //glRectd(0,0,1920,1080);
    int cw=15,l=0;
    for(int j=0;l<coord.size();j++){
      glBegin(GL_LINES);
      for(int i=0;i<cw;i++){
	glVertex2d(coord[i+cw*j].x1,
		   coord[i+cw*j].y1);
	l++;
      }
      glEnd();
    }
    glPopMatrix();
    
    SDL_GL_SwapBuffers();
    SDL_Delay(64);
    
    
    SDL_Event e;							
    while(SDL_PollEvent(&e))		
      if(e.type==SDL_KEYDOWN && e.key.keysym.sym==SDLK_ESCAPE)	
	goto end;						
  }	
 end:
  return 0;
}
