/*!  \file writeData.cpp
 *   \brief Write out the solution
 *
 *   \author Julian Pittard (Original version 09.08.11)
 *   \version 0.1-development (Undómiel):
 *   \date Last modified: 09.08.11 (JMP)
 */

#undef MAIN

//  Header files:
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "constants.h"
#include "gas.h"
#include "global.h"
#include "grid.h"
#include "physconst.h"
#include "star.h"

#include "fitsio.h"      //  For FITSIO (from cfitsio directory)
#include "longnam.h"     //  For FITSIO (from cfitsio directory)
// File scope
int bitpix = FLOAT_IMG;                 //  Single precision.

// Function declarations
string generate_asciifilename(int nfile);
string generate_fitsfilename(int nfile);
void FitsPrintError( int status );
void writeFitsKeywords(fitsfile *fptr);
void writeFitsKeywordsCWB(fitsfile *fitsptr);
void write2dFitsImage(fitsfile *fptr, char* extname, long naxis, long naxes[2], double array2d[][imax]);
void write2dFitsImage(fitsfile *fptr, char* extname, long naxis, long naxes[2], float array1d[]);
void write2dFitsImage(fitsfile *fptr, char* extname, long nzones, double array1d[]);
void write2dFitsImage(fitsfile *fptr, char* extname, long nzones, float array1d[]);
extern void quit();

using namespace std;


void prin1D(){

  string outfile;
  int i=0;
  int j=0;
  int k=0;
  int m=0;

  // Write out data 
  outfile = generate_asciifilename(nfile1d);// generate filename
  //cout << outfile << "\n";
  ofstream file(outfile.c_str()); //C++ method (outfile must be a string type)
  if (!file) {
    cerr << "Failed to open output file: " << outfile;
    quit();
  } 
  file << "# Problem = " << problem << "; ncycle = " << ncycle << "; t = " << t << "\n";
  file << "# zxa             rho             pre              vx             col\n";
  file << scientific << setprecision(4);
  for (int i = lg.irs; i <= lg.ire; ++i){
    file << lg.zxa[i] << "\t" << lg.P0[iqd][k][j][i] << "\t" << lg.P0[iqe][k][j][i] << "\t" << lg.P0[iqu0][k][j][i] << "\t" << lg.P0[iqal0][k][j][i] << "\n";
  }
  file.close();
  ++nfile1d;
  return;
}


void prin2D(){

  // Write out a file containing x,y,rho,pressure,vx,vy.
  // To plot in gnuplot use:
  // set term png
  // set output "sod_0000_2d_0009.png"
  // set size square
  // set logscale cb
  // set cbrange [0.1:15.0]
  // plot 'sod_0000_2d_0009.ascii' us 1:2:4 with image pixels
  // set term x11
  
  string outfile;
  int i=0;
  int j=0;
  int k=0;
  int m=0;

  // Write out data 
  outfile = generate_asciifilename(nfile2d);// generate filename
  //cout << outfile << "\n";
  ofstream file(outfile.c_str()); //C++ method (outfile must be a string type)
  if (!file) {
    cerr << "Failed to open output file: " << outfile;
    quit();
  } 
  file << "# Problem = " << problem << "; ncycle = " << ncycle << "; t = " << t << "\n";
  file << "# zxa             zya              rho             pre              vx             col\n";
  file << scientific << setprecision(4);
  for (int j = lg.jrs; j <= lg.jre; j++){
    for (int i = lg.irs; i <= lg.ire; ++i){
      file << lg.zxc[i] << "\t" << lg.zyc[j] << "\t" << lg.P0[iqd][k][j][i] << "\t" << lg.P0[iqe][k][j][i] << "\t" << lg.P0[iqu0][k][j][i] << "\t" << lg.P0[iqal0][k][j][i] << "\n";
    }
    file << "\n";
  }
  file.close();
  ++nfile2d;
  return;
}


void prin2Dfits(){

  // Write out a 2D FITS file compatible with the X-ray code
  
  int i,ndim;
  char type[20];
  char cnum[5];             // needs an extra space for end of string character

  float scrh;
  FILE *fptr;
  
  //  FITSIO variables.
  int simple = 1;
  int extend = 1;
  int status = 0;                         //  FITSIO status.
  fitsfile *fitsptr;                      //  Pointer to fits file.
  int bitpix = FLOAT_IMG;                 //  Single precision.
  long naxis = 2;
  long naxes[2];
  long group = 1;
  long fpixel[2];
  long nelements;
  char extname[FLEN_VALUE];
  
  string outfile;    
  const char *filename;     //  Character to store filename

  const char* ptrbuf1;      //  To call FITS routines expecting a const char
  char buf1[30];            //  without an error

  double data2d[jmax][imax];

  ndim = 2;
  outfile = generate_fitsfilename(nfile2d);     // generate filename
  filename = outfile.c_str();                   // convert to character array

  // Delete any existing file of that name and write out to file
  remove( filename );
  if ( fits_create_file(&fitsptr, filename, &status) ){ 
    FitsPrintError( status );
  }

  naxes[0] = imax;
  naxes[1] = jmax;

  //  Create the density image.
  if ( fits_create_img(fitsptr, bitpix, naxis, naxes, &status) ){
    cerr << "  Error creating rho data image" << endl;
    FitsPrintError( status );
  }

  //  Write keywords
  writeFitsKeywords(fitsptr);

  // Write the first (density) image
  //  Put the data into a 1-D array  
  float *ptr_flt = new float[imax*jmax];
  int n = 0;
  for (int j=0; j < jmax; j++){
    for (int i=0; i < imax; i++){
      ptr_flt[n] = lg.P0[iqd][0][j+nghost][i+nghost];
      //cout << ptr_flt[i] << "\t";
      n++;
    }
  }
    
  //  Write the 1-D array, cunningly stored in the same way as a 2-D
  //  image, into the FITS image.
  if ( fits_write_2d_flt(fitsptr, group, naxes[0], 
                           naxes[0], naxes[1], ptr_flt, &status) ){ 
    FitsPrintError( status );
  } 
  
  delete []ptr_flt; // clean up array

  // Create and write the other data as (2d) images
  for (int j=0; j < jmax; j++){
    for (int i=0; i < imax; i++){
      data2d[j][i] = lg.P0[iqe][0][j+nghost][i+nghost];
    }
  }
  write2dFitsImage(fitsptr,"pressure",naxis,naxes,data2d);

  for (int j=0; j < jmax; j++){
    for (int i=0; i < imax; i++){
      data2d[j][i] = lg.P0[iqu0][0][j+nghost][i+nghost];
    }
  }
  write2dFitsImage(fitsptr,"vx",naxis,naxes,ptr_flt);

  for (int j=0; j < jmax; j++){
    for (int i=0; i < imax; i++){
      data2d[j][i] = lg.P0[iqu0+1][0][j+nghost][i+nghost];
    }
  }
  write2dFitsImage(fitsptr,"vy",naxis,naxes,ptr_flt);

  //write2dFitsImage(fitsptr,"vz",naxis,naxes,zuz);
  
  for (int m = 0; m < smax; m++){
    sprintf(cnum, "%4d",m);
    for (int i = 0; i < 4; i++){
      if (cnum[i] == ' ') cnum[i] = '0';
    }
    //strcpy(type,"as");
    //strcat(type,cnum);
    strcpy(type,"colour");
    for (int j=0; j < jmax; j++){
      for (int i=0; i < imax; i++){
        data2d[j][i] = lg.P0[iqal0+m][0][j+nghost][i+nghost];
      }
    }
    write2dFitsImage(fitsptr,type,naxis,naxes,ptr_flt); //advected scalars
  }
    
  //Write axis data as images
  write2dFitsImage(fitsptr,"zxa",naxes[0]+1,lg.zxa);
  write2dFitsImage(fitsptr,"zya",naxes[1]+1,lg.zya);
  write2dFitsImage(fitsptr,"zdx",naxes[0]+1,lg.zdx);
  write2dFitsImage(fitsptr,"zdy",naxes[1]+1,lg.zdy);

  //  Close the output fits file.
  if ( fits_close_file(fitsptr, &status) ){   //  Close the file.
    FitsPrintError( status ); 
  }
  //cout << " Finished writing " << outfile << endl;

  // Check for any error, and if so print out error messages
  if (status > 0){
    FitsPrintError(status);
    quit();
  }

  ++nfile2d;
  return;
}


void FitsPrintError( int status ){
  if (status)
  {
    fits_report_error(stderr, status);    //  Print error report.
    //exit( status );    //  Terminate the program, returning error status.
    quit();    //  Terminate the program
  }
  return;

}


void writeFitsKeywords(fitsfile *fitsptr){

  //Locals
  int ndim = nd;
  int c_max = cmax;
  int s_max = smax;
  double d_floor = dfloor;
  double p_floor = pfloor;

  int status = 0;                         //  FITSIO status.

  char extname[FLEN_VALUE];
  double gamma,courant,avg_mass;

  char pre_fix[30];
  char pre_fixstats[70];
  char chydrosolver[70]; //  character array to store hydro solver
  char cproblem[70];

  //  Write keywords, in particular the coordinate keywords.
  //  First write current date and time (UTC).
  if ( fits_write_date(fitsptr, &status ) ){
    cerr << "  Error writing DATE keyword."  << endl;
    FitsPrintError( status );
  }

  //  Now write the density extension name
  strcpy(extname,"density");
  if ( fits_update_key(fitsptr, TSTRING, "VH_DATA", extname, 
			 "Extension name", &status ) ){
    cerr << "  Error writing EXTNAME version keyword."  << endl;
    FitsPrintError( status );
  }

  //Write introductory (mostly logical) keywords
  strcpy(cproblem,problem.c_str());   // copy string into char array
  fits_update_key(fitsptr,TSTRING,"problem",cproblem,"problem",&status);

  //Write integers...
  fits_update_key(fitsptr,TINT,"ncycle",&ncycle,"ncycle",&status);

  // Now the doubles. const parameters must first be copied into non-const variables
  // before they can be output
  float avgmass1 = 1.0e-24;
  float avgmass2 = 1.0e-24;
  
  fits_update_key(fitsptr,TFLOAT,"avgmass1",&avgmass1,"average mass 1",&status);
  fits_update_key(fitsptr,TFLOAT,"avgmass2",&avgmass2,"average mass 2",&status);

  // Problem specific code to come here...
  //

  if (problem == "CWB"){
    writeFitsKeywordsCWB(fitsptr);
  }

  // Check status before exiting

  if (status != 0){
    cerr << "  Error writing CTYPE keywords."  << endl;
    FitsPrintError( status );
  }

  //cout << "writeFitsKeywords: Keywords written successfully\n";
  return;
}


void writeFitsKeywordsCWB(fitsfile *fitsptr){

  int status = 0;
  int intval;
  
  fits_update_key(fitsptr,TDOUBLE,"xpos_1",&stars[0].xpos,"(cm)",&status);
  fits_update_key(fitsptr,TDOUBLE,"ypos_1",&stars[0].ypos,"(cm)",&status);
  //fits_update_key(fitsptr,TDOUBLE,"zpos1",&stzpos[0],"(cm)",&status);
  fits_update_key(fitsptr,TDOUBLE,"xpos_2",&stars[1].xpos,"(cm)",&status);
  fits_update_key(fitsptr,TDOUBLE,"ypos_2",&stars[1].ypos,"(cm)",&status);

  double X1 = 0.705;
  double X2 = 0.705;
  double Y1 = 0.275;
  double Y2 = 0.275;
  
  fits_update_key(fitsptr,TDOUBLE,"X1",&X1,"X1",&status);
  fits_update_key(fitsptr,TDOUBLE,"X2",&X2,"X2",&status);
  fits_update_key(fitsptr,TDOUBLE,"Y1",&Y1,"Y1",&status);
  fits_update_key(fitsptr,TDOUBLE,"Y2",&Y2,"Y2",&status);

  // Check status before exiting

  if (status != 0){
    cerr << "  Error writing CWB keywords."  << endl;
    FitsPrintError( status );
  }

  //cout << "writeFitsKeywordsCWB: Keywords written successfully\n";
  return;
}


void write2dFitsImage(fitsfile *fitsptr, char* extname, long naxis, long naxes[2], double array2d[][imax]){

  int status = 0;
  long group = 1;

  //  Create the new image.
  if ( fits_create_img(fitsptr, bitpix, naxis, naxes, &status) ){
    cerr << "  Error creating " << extname << " data image" << endl;
    FitsPrintError( status );
  }
  if ( fits_update_key(fitsptr, TSTRING,"VH_DATA", extname, 
			 "Extension name", &status ) ){
    cerr << "  Error writing EXTNAME version keyword."  << endl;
    FitsPrintError( status );
  }
     
  //  Get the passed data into a 1-D array 
  // NOTE: Not yet coded for 2d slices not in the xy plane 
  float *ptr_flt = new float[imax*jmax];
  int n = 0;
  for (int j=0; j < jmax; j++){
    for (int i=0; i < imax; i++){
      ptr_flt[n] = array2d[j][i];
      n++;
    }
  }
    
  //  Write the 1-D array, cunningly stored in the same way as a 2-D
  //  image, into the FITS image.
  if ( fits_write_2d_flt(fitsptr, group, naxes[0], 
			   naxes[0], naxes[1], ptr_flt, &status) ){ 
    FitsPrintError( status );
  } 

  // Write pixel keywords
  //writeFitsKeywordsPixels(fitsptr);

  delete []ptr_flt; // clean up array

  return;
}


void write2dFitsImage(fitsfile *fitsptr, char* extname, long naxis, long naxes[2], float array1d[]){

  int status = 0;
  long group = 1;

  //  Create the new image.
  if ( fits_create_img(fitsptr, bitpix, naxis, naxes, &status) ){
    cerr << "  Error creating " << extname << " data image" << endl;
    FitsPrintError( status );
  }
  if ( fits_update_key(fitsptr, TSTRING,"VH_DATA", extname, 
			 "Extension name", &status ) ){
    cerr << "  Error writing EXTNAME version keyword."  << endl;
    FitsPrintError( status );
  }
     
  //  Get the passed data into a 1-D array  
  float *ptr_flt = new float[naxes[0]];
  for (int i=0; i < naxes[0]; i++){
    ptr_flt[i] = array1d[i];
  }
    
  //  Write the 1-D array, cunningly stored in the same way as a 2-D
  //  image, into the FITS image.
  if ( fits_write_2d_flt(fitsptr, group, naxes[0], 
			   naxes[0], naxes[1], ptr_flt, &status) ){ 
    FitsPrintError( status );
  } 

  delete []ptr_flt; // clean up array

  return;
}


void write2dFitsImage(fitsfile *fitsptr, char* extname, long nzones, float array1d[]){

  int status = 0;
  long group = 1;
  long naxis;
  long naxes[2];

  naxis = 2;
  naxes[0] = nzones;
  naxes[1] = 1;

  //  Create the new image.
  if ( fits_create_img(fitsptr, bitpix, naxis, naxes, &status) ){
    cerr << "  Error creating " << extname << " data image" << endl;
    FitsPrintError( status );
  }
  if ( fits_update_key(fitsptr, TSTRING,"VH_DATA", extname, 
			 "Extension name", &status ) ){
    cerr << "  Error writing EXTNAME version keyword."  << endl;
    FitsPrintError( status );
  }
     
  //  Get the passed data into a 1-D array  
  float *ptr_flt = new float[nzones];
  for (int i=0; i < nzones; i++){
    ptr_flt[i] = array1d[i];
  }
    
  //  Write the 1-D array, cunningly stored in the same way as a 2-D
  //  image, into the FITS image.
  if ( fits_write_2d_flt(fitsptr, group, naxes[0], 
			   naxes[0], naxes[1], ptr_flt, &status) ){ 
    FitsPrintError( status );
  } 

  delete []ptr_flt; // clean up array

  return;
}


void write2dFitsImage(fitsfile *fitsptr, char* extname, long nzones, double array1d[]){

  int status = 0;
  long group = 1;
  long naxis;
  long naxes[2];

  naxis = 2;
  naxes[0] = nzones;
  naxes[1] = 1;

  //  Create the new image.
  if ( fits_create_img(fitsptr, bitpix, naxis, naxes, &status) ){
    cerr << "  Error creating " << extname << " data image" << endl;
    FitsPrintError( status );
  }
  if ( fits_update_key(fitsptr, TSTRING,"VH_DATA", extname, 
			 "Extension name", &status ) ){
    cerr << "  Error writing EXTNAME version keyword."  << endl;
    FitsPrintError( status );
  }
     
  //  Get the passed data into a 1-D array  
  float *ptr_flt = new float[nzones];
  for (int i=0; i < nzones; i++){
    ptr_flt[i] = array1d[i];
  }
    
  //  Write the 1-D array, cunningly stored in the same way as a 2-D
  //  image, into the FITS image.
  if ( fits_write_2d_flt(fitsptr, group, naxes[0], 
			   naxes[0], naxes[1], ptr_flt, &status) ){ 
    FitsPrintError( status );
  } 

  delete []ptr_flt; // clean up array

  return;
}


string generate_asciifilename(int nfile)
//Generate the name for ASCII file dumps. Always include the processor number
//(for serial code this is 0000).
{
  char cnum[5];             // needs an extra space for end of string character
  char ndnum[2];            // needs an extra space for end of string character
  string srank,snfile,sdim;
  string filename;

  // Create strings for procRank and nfile
  sprintf(cnum, "%4d",procRank);
  for (int i = 0; i < 4; i++){
    if (cnum[i] == ' ') cnum[i] = '0';
  }
  srank = cnum;

  sprintf(cnum, "%4d",nfile);
  for (int i = 0; i < 4; i++){
    if (cnum[i] == ' ') cnum[i] = '0';
  }
  snfile = cnum;

  //sdim = to_string(nd);
  sprintf(ndnum, "%1d",nd);
  sdim = ndnum;
  sdim += 'd';
 
  //The ASCII filename takes the form: prefx_AAAA_BB_CCCC.ascii
  // where AAAA = rank number of processor (e.g. 0000)
  //         BB = dimensionality of data   (e.g. 1d)
  //       CCCC = file number in sequence  (e.g. 0010)
  
  filename = prefx + '_' + srank + '_' + sdim + '_' + snfile + ".ascii";

  return filename;
}


string generate_fitsfilename(int nfile)
//Generate the name for FITS file dumps. Always include the processor number
//(for serial code this is 0000).
{
  char cnum[5];             // needs an extra space for end of string character
  char ndnum[2];            // needs an extra space for end of string character
  string srank,snfile,sdim;
  string filename;

  // Create strings for procRank and nfile
  sprintf(cnum, "%4d",procRank);
  for (int i = 0; i < 4; i++){
    if (cnum[i] == ' ') cnum[i] = '0';
  }
  srank = cnum;

  sprintf(cnum, "%4d",nfile);
  for (int i = 0; i < 4; i++){
    if (cnum[i] == ' ') cnum[i] = '0';
  }
  snfile = cnum;

  //sdim = to_string(nd);
  sprintf(ndnum, "%1d",nd);
  sdim = ndnum;
  sdim += 'd';
 
  //The ASCII filename takes the form: prefx_AAAA_BB_CCCC.ascii
  // where AAAA = rank number of processor (e.g. 0000)
  //         BB = dimensionality of data   (e.g. 1d)
  //       CCCC = file number in sequence  (e.g. 0010)
  
  filename = prefx + '_' + srank + '_' + sdim + '_' + snfile + ".fits";

  return filename;
}
