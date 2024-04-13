/*********************************************************************
  Lucas and Kanade, 1981
  With non-hierarchical modifications using Simoncelli, 
    Adelson and Heeger, CVPR 91
  cc -g lucas.c -lm -o lucas
**********************************************************************/

float         floatpic[FIVE][PIC_X][PIC_Y];
float         Ix[PIC_X][PIC_Y],
              Iy[PIC_X][PIC_Y],
              It[PIC_X][PIC_Y],
              Imag[PIC_X][PIC_Y];
float         full_vels[2][PIC_X][PIC_Y];
int           int_size_x,
              int_size_y;
int           pic_x,
              pic_y,
              pic_t;
int           binary,
              raw_statistics,
              output_smooth;
float         raw_mag,
              MAX_raw_mag;
int           fd_temp[8];
float         actual_x,
              actual_y,
              size_x,
              size_y,
              offset_x,
              offset_y;
float         diff1[2],
              diff2[2];
float         data[8][150][150];
int           check_eigen_calc();
float         diff_x(),diff_y(),diff_t();
float         PsiER(),PsiEN(),norm();


/*********************************************************************
 * Usage
 *********************************************************************/

void usage(progname)
char *progname;
{
   printf(USAGE, progname);
   printf("Read enough images so that smoothing with a 3D Gaussian\n");
   printf(" with <sigma> leaves 5 smoothed images centered at ");
   printf("<central file number>\n");
   printf("sigma is the standard deviation used in smoothing the images\n");
   printf("sigma==0.0 means no smoothing: 5 unsmoothed images are used\n");
   printf("<input path> - directory where input data resides\n");
   printf("<output path> - directory where computed flow fields are put\n");
   printf("-L - standard Lucas and Kanade with presmoothing and thresholding\n");
   printf("-M - modified Lucas and Kanade (Simoncelli, Adelson and Heeger)\n");
   printf("default -L if any flag other than -M specified or no flag specified\n");
   printf("-B <cols> <rows> - use binary instead of black and white rasterfiles as image input\n");
   printf("-S <smoothed path> - creates appropriate files for the smoothed input data\n");
   printf("   and writes to the specifed directory\n");
   printf("-C <correct filename> - perform error analysis using correct velocity field\n");
   printf("-T <float> - specify a magnitude threshold for spatial derivatives\n");
   printf("             for the raw normal velocity data\n");
   printf("If -T 0.0 is specified histogram and cumulative statistics are printed\n");
   printf("If -T is not specified default value is 0.0 with no statistics printed\n");
}

/*********************************************************************/
/* Read input images                                                 */
/*********************************************************************/

/* Function to compose the file name with a frame number and extension */


char *ComposeFname(s, nr)
char *s;
int nr;
{
   char r[4], fnam[64];
   int l;

   strcpy(r,"000");
   for (l=(strlen(r)-1);l;l--)
   {
      r[l] += (char) nr % 10;
      nr /= 10;
   }
   sprintf(fnam,"%s%s%s",s,r,".raw");
   strcpy(s, fnam);
   return(s);
}

/* Function to read one file with name and frame number */

void readfile (fname, inpic, pic_x, pic_y, header)
char fname[32];
unsigned char inpic[PIC_X][PIC_Y];
int *pic_x, *pic_y;
unsigned char header[HEAD];
{
   static int once=TRUE;
   int j, fd, no_bytes = 0;
   int ints[8];

   if ((fd=open (fname,O_RDONLY)) >0)
   {
      if (!binary)
      {
         if (once)
         {
            no_bytes += read (fd, ints, HEAD);
            (*pic_y) = ints[1];
            (*pic_x) = ints[2];
            if((*pic_x) > PIC_X || (*pic_y) > PIC_Y)
            {
               printf("Fatal error: images are too big.\n");
               printf("Use only images with %d columns, %d rows\n",
                      PIC_X, PIC_Y);
               exit (1);
            }
            once = FALSE;
         }
         else no_bytes += read (fd, header, HEAD);
      }
      for (j=0; j<(*pic_x); j++)
         no_bytes += read (fd, &inpic[j][0], (*pic_y));
      printf ("File %s read (%d bytes)\n", fname, no_bytes);
      close (fd);
   }
   else 
   {
      printf ("File %s does not exist in readfiles.\n", fname);
      exit (1);
   }
}

void readfiles(path,stem,pic,sigma,pic_t,pic_x,pic_y,start,end,header)
char path[32], stem[32];
unsigned char pic[PIC_T][PIC_X][PIC_Y];
float sigma;
int pic_t, *pic_x, *pic_y, start, end;
unsigned char header[HEAD];
{
   char fname[32];
   int i, time;

   strcpy (fname, "");
   printf ("Reading Files...\n");
   if ((end-start) < 0)
   {
      printf ("\nSpecified time for writing file incorrect\n");
      exit (1);
   }
   for (i=start, time=0; i<=end; i++, time++)
   {
      sprintf (fname, "%s/%s", path, stem);
      ComposeFname (fname, i);
      readfile (fname, pic[time], pic_x, pic_y, header);
   }
   fflush (stdout);
} /* End of readfiles */


/*********************************************************************/
/*   Write smoothed images                                           */
/*********************************************************************/
void writefiles(path,stem,result,sigma,pic_t,pic_x,pic_y,start,end,header)
char path[32], stem[32];
float result[FIVE][PIC_X][PIC_Y];
float sigma;
int pic_t, pic_x, pic_y, start, end;
unsigned char header[HEAD];
{
   char fname[32];
   int i, j, fd, time, no_bytes;
   unsigned char pic[PIC_X][PIC_Y];
   float term;

   if(sigma==0.0) return;
   printf("\nWriting smoothed files...\n");

   for (time=start; time<=end; time++)
   {
      for (i=0; i<pic_x; i++)
         for (j=0; j<pic_y; j++)
         {
            term = result[time-start][i][j];
            if(term > 255.0) term = 255.0;
            if(term < 0.0) term = 0.0;
            pic[i][j] = (int) (term+0.5);
         }
      no_bytes = 0;
      sprintf (fname, "%s/smoothed.%s%d-%3.1f", path, stem, time, sigma);
      if ((fd=creat (fname,0644))!=(-1))
      {
         no_bytes += write (fd,&header[0],HEAD);/*Write 32 byte raster header*/
         for (j=0; j<pic_x; j++)
            no_bytes += write (fd, &pic[j][0], pic_y);
         printf ("File %s written (%d bytes)\n", fname, no_bytes);
      }
      else 
      {
         printf ("File %s cannot be written\n", fname);
         exit (1);
      }
      fflush (stdout);
      close (fd);
   }
} /* End of writefiles */


/*********************************************************************/
/* Smooth the image sequence using a 3D separable Gaussian filter.   */
/* Perform 3D Gaussian smoothing by separable convolution            */
/*********************************************************************/
void Smooth3D_Gaussian (path, stem, floatpic, sigma, pic_t, pic_x, pic_y,
                   start, middle, end, header)
char path[32], stem[32];
float floatpic[FIVE][PIC_X][PIC_Y];
float sigma;
int pic_x, pic_y, pic_t, start, middle, end;
unsigned char header[HEAD];
{
   char fname[100];
   float mask[100], term, sum;
   int time, size, first, last, i, j, k, offset, a, b;
   unsigned char inpic[PIC_X][PIC_Y];
   float pic[PIC_X][PIC_Y];

   pic_t = end-start+1;
   size = (int) 6*sigma+1;
   if (size%2==0) size = size+1;
   offset = size/2;
   if (pic_t < size)
   {
      printf("\nFatal error: not enough images\n");
      exit(1);
   }
   if (sigma != 0.0)
   {
      sum = 0.0;
      for (i=0; i<size; i++)
      {
         mask[i] = (1.0/(sqrt(2.0*3.141592654)*sigma))*
             exp(-(i-offset)*(i-offset)/(2.0*sigma*sigma));
         sum += mask[i];
      }
   }
   else { mask[0] = 1.0; sum = 1.0; }

   if (output_smooth)
   {
      printf("Start: %d End: %d Middle: %d\n",start,end,middle);
      printf("\nStart of 3D convolution\n");
      printf("Size: %d Offset: %d\nMask values: ",size,offset);
      for (i=0; i<size; i++)
         printf ("%f ",mask[i]);
      printf ("\nSum of mask values: %f\n",sum);
   }
   if (sigma != 0.0)
   {
      for (i=0;i<pic_x;i++)
         for (j=0;j<pic_y;j++)
         {
            pic[i][j] = 0.0;
            for (k=0; k<FIVE; k++)
               floatpic[k][i][j] = 0.0;
         }
      for (time = start; time <= end; time++)
      {
         sprintf (fname, "%s/%s", path, stem);
         ComposeFname (fname, time);
         readfile (fname, inpic, &pic_x, &pic_y, header);

         first = time-start < size ? 0 : time-start-size+1;
         last  = time-start < pic_t-size? time-start : pic_t-size;
         for (k = first; k <= last ; k++)
         {
            if (output_smooth) printf ("Time: %d Frame: %d\n", time, k);
            for(i=0;i<pic_x;i++)
               for(j=0;j<pic_y;j++)
               {
                  floatpic[k][i][j] += (inpic[i][j]*mask[time-start-k]);
               }
         }
      }
      if(output_smooth) printf("Convolution in t direction completed\n");
      for (k = 0; k < FIVE; k++)
      {
         if (output_smooth) printf ("Frame: %d\n", k+offset);

         for(i=0;i<offset;i++)
            for(j=0;j<pic_y;j++)
               pic[i][j] = floatpic[k][i][j];

         for(i=offset;i<pic_x-offset;i++) {
            for(j=0;j<offset;j++)
               pic[i][j] = floatpic[k][i][j];
            for(j=offset;j<pic_y-offset;j++)
            {
               term = 0.0;
               for(a=-offset;a<=offset;a++)
               {
                  term = term + (floatpic[k][i+a][j]*mask[a+offset]);

               }
               pic[i][j] = term;
            }
            for(j=pic_y-offset;j<pic_y;j++)
               pic[i][j] = floatpic[k][i][j];
	 }

         for(i=pic_x-offset;i<pic_x;i++)
            for(j=0;j<pic_y;j++)
               pic[i][j] = floatpic[k][i][j];

         if(output_smooth) printf("Convolution in x direction completed\n");
   
         for(i=offset;i<pic_x-offset;i++)
            for(j=offset;j<pic_y-offset;j++)
            {
               term = 0.0;
               for(b=-offset;b<=offset;b++)
               {
                  term = term + (pic[i][j+b])*mask[b+offset];
               }
               floatpic[k][i][j] = term;
            }
         if(output_smooth) printf("Convolution in y direction completed\n");
      }
      if(output_smooth) printf("End of Convolution\n");
      printf("Input images smoothed\n");
   }
   else /* No smoothing */
   {
      for (time = start; time <= end; time++)
      {
         sprintf (fname, "%s/%s", path, stem);
         ComposeFname (fname, time);
         readfile (fname, inpic, &pic_x, &pic_y, header);

         for(i=0;i<pic_x;i++)
            for(j=0;j<pic_y;j++)
               floatpic[time-start][i][j] = (float) inpic[i][j];
      }
   }
   fflush(stdout);
}


/*********************************************************************/
/* Compute spatio-temporal derivatives by applying a 4 point         */
/* central difference kernel to each direction x, y and t.           */
/*********************************************************************/

/***********************************************************
   Compute a 4 point central difference kernel             
***********************************************************/
void calc_diff_kernel(diff_kernel)
float diff_kernel[5];
{
   diff_kernel[0] = -1.0/12.0;
   diff_kernel[1] = 8.0/12.0;
   diff_kernel[2] = 0.0;
   diff_kernel[3] = -8.0/12.0;
   diff_kernel[4] = 1.0/12.0;
}


/************************************************************
   Apply 1D real kernels in the x direction
************************************************************/
float diff_x(floatpic,kernel,x,y,n)
float floatpic[FIVE][PIC_X][PIC_Y];
float kernel[5];
int x,y,n;
{
   int i;
   float sum;

   sum = 0.0;
   for(i=(-n);i<=n;i++)
   {
      sum += kernel[i+2]*floatpic[2][x+i][y];
   }
   return(sum);
}

/************************************************************
   Apply 1D real kernels in the y direction
************************************************************/
float diff_y(floatpic,kernel,x,y,n)
float floatpic[FIVE][PIC_X][PIC_Y];
float kernel[5];
int x,y,n;
{
   int i;
   float sum;

   sum = 0.0;
   for(i=(-n);i<=n;i++)
   {
      sum += kernel[i+2]*floatpic[2][x][y+i];
   }
   return(sum);
}


/************************************************************
   Apply 1D real kernels in the t direction.
************************************************************/
float diff_t(floatpic,kernel,x,y,n)
float floatpic[FIVE][PIC_X][PIC_Y];
float kernel[5];
int x,y,n;
{
   int i;
   float sum;

   sum = 0.0;
   for(i=(-n);i<=n;i++)
   {
      sum += kernel[i+2]*floatpic[2+i][x][y];
   }
   return(sum);
}


/************************************************************
   Compute spatio-temporal derivatives                     
************************************************************/
void compute_ders (floatpic, Ix, Iy, It, pic_t, pic_x, pic_y, n)
float floatpic[FIVE][PIC_X][PIC_Y];
float Ix[PIC_X][PIC_Y];
float Iy[PIC_X][PIC_Y];
float It[PIC_X][PIC_Y];
int pic_t, pic_x, pic_y, n;
{
   int i,j;
   float kernel[5];

   for(i=0;i<PIC_X;i++)
      for(j=0;j<PIC_X;j++)
      {
         Ix[i][j] = Iy[i][j] = It[i][j] = 0.0;
      }
   calc_diff_kernel(kernel);
   for(i=n;i<pic_x-n;i++)
      for(j=n;j<pic_y-n;j++)
      {
         Ix[i][j] = diff_x(floatpic,kernel,i,j,2);
         Iy[i][j] = diff_y(floatpic,kernel,i,j,2);
         It[i][j] = diff_t(floatpic,kernel,i,j,2);
      }
   printf("Spatio-Temporal Intensity Derivatives Computed\n");
   fflush(stdout);
}


/******************************************************************/
/*   Compute velocities                                           */
/******************************************************************/

/******************************************************************/
/* Compute all eigenvalues and eigenvectors of a real symmetric   */
/* matrix a[N][N]. On output elements of a above the diagonal are */
/* destroyed. d[N] returns the eigenvalues of a. v[N][N] is a     */
/* matrix whose columns contain, on output, the normalized        */
/* eigenvectors of a. nrot returns the number of Jacobi rotations */
/* that were required.                                            */
/******************************************************************/

/**********************************************************
   Do rotations required by Jacobi Transformation         
**********************************************************/
void rotate(a,i,j,k,l,h,g,s,tau)
float a[N][N],s,tau;
int i,j,k,l;
float *h,*g;
{
   (*g) = a[i][j];
   (*h) = a[k][l];
   a[i][j] = (*g)-s*((*h)+(*g)*tau);
   a[k][l] = (*h)+s*((*g)-(*h)*tau);
}

/*********************************************************
   Compute Jacobi transformation
*********************************************************/
void jacobi(aa,n,d,v,nrot)
float aa[N][N],d[N],v[N][N];
int n,*nrot;
{
   int j,iq,ip,i;
   float thresh,theta,tau,t,sm,s,h,g,c;
   float b[N],z[N],a[N][N];

   if(n!=N) {
      fprintf(stderr,"\nFatal error: n not N in jacobi\n",N);
      exit(1);
   }
   for(ip=0;ip<n;ip++) /* Initialize to the identity matrix */
   {
      for(iq=0;iq<n;iq++) v[ip][iq] = 0.0;
      for(iq=0;iq<n;iq++) a[ip][iq] = aa[ip][iq]; /* Don't destroy aa */
      v[ip][ip] = 1.0;
   }
   /* Initialize b and d to the diagonals of a */
   for(ip=0;ip<n;ip++)
   {
      b[ip] = d[ip] = a[ip][ip];
      z[ip] = 0.0;
   }
   *nrot = 0;
   for(i=0;i<100;i++)
   {
      sm = 0.0;
      for(ip=0;ip<(n-1);ip++)
      {
         for(iq=ip+1;iq<n;iq++)
            sm += fabs(a[ip][iq]);
      }

      /* Normal return, which relies on quadratic convergence to
     machine underflow */
      if(sm == 0.0) return;

      if(i<3) thresh=0.2*sm/(n*n); /* on the first three sweeps */
      else thresh = 0.0; /* the rest of the sweeps */

      for(ip=0;ip<(n-1);ip++)
      {
         for(iq=ip+1;iq<n;iq++)
         {
            g = 100.0*fabs(a[ip][iq]);
            /* After 4 sweeps skip the rotation if the
         off diagonal element is small */
            if(i>3 && fabs(d[ip])+g == fabs(d[ip])
                && fabs(d[iq])+g == fabs(d[iq])) a[ip][iq] = 0.0;
            else if(fabs(a[ip][iq]) > thresh)
            {
               h = d[iq]-d[ip];
               if(fabs(h)+g==fabs(h)) t=(a[ip][iq])/h;
               else 
               {
                  theta = 0.5*h/(a[ip][iq]);
                  t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                  if(theta < 0.0) t = -t;
               }
               c = 1.0/sqrt(1.0+t*t);
               s = t*c;
               tau = s/(1.0+c);
               h = t*a[ip][iq];
               z[ip] -= h;
               z[iq] += h;
               d[ip] -= h;
               d[iq] += h;
               a[ip][iq] = 0.0;
               for(j=0;j<ip-1;j++)
                  rotate(a,j,ip,j,iq,&h,&g,s,tau);
               for(j=ip+1;j<iq-1;j++)
                  rotate(a,ip,j,j,iq,&h,&g,s,tau);
               for(j=iq+1;j<n;j++)
                  rotate(a,ip,j,iq,j,&h,&g,s,tau);
               for(j=0;j<n;j++)
                  rotate(v,j,ip,j,iq,&h,&g,s,tau);
               ++(*nrot);
            }
         }
      }
      for(ip=0;ip<n;ip++)
      {
         b[ip] += z[ip];
         d[ip] = b[ip];
         z[ip] = 0.0;
      }
   }
   /* fprintf(stderr,"\nFatal error: too many iterations in jacobi\n"); */
}

/*******************************************************************
   Check eigenvector and eigenvalue computation for 2*2 matrix     
*******************************************************************/
int check_eigen_calc(mm,d,v,n,diff1,diff2,length1,length2,angle)
float mm[N][N],d[N],v[N][N],diff1[N],diff2[N],*angle,*length1,*length2;
int n;
{
   int status;

   status = TRUE;
   /* Compute angle between two eigenvectors - should be orthogonal */
   (*angle)=acos((v[0][0]*v[0][1]+v[1][0]*v[1][1])/
       (sqrt(v[0][0]*v[0][0]+v[1][0]*v[1][0])*
       sqrt(v[0][1]*v[0][1]+v[1][1]*v[1][1])))*180.0/3.1415926535;
   if((*angle) < 89.5 && (*angle) > 90.5)
   {
      status = FALSE;
   }

   /* Eigenvector test */
   diff1[0] = mm[0][0]*v[0][0]+mm[0][1]*v[1][0];
   diff1[1] = mm[1][0]*v[0][0]+mm[1][1]*v[1][0];
   diff1[0] = diff1[0] - d[0]*v[0][0];
   diff1[1] = diff1[1] - d[0]*v[1][0];
   if(((*length1)=sqrt(diff1[0]*diff1[0]+diff1[1]*diff1[1])) > 0.1)
   {
      status = FALSE;
   }
   diff2[0] = mm[0][0]*v[0][1]+mm[0][1]*v[1][1];
   diff2[1] = mm[1][0]*v[0][1]+mm[1][1]*v[1][1];
   diff2[0] = diff2[0] - d[1]*v[0][1];
   diff2[1] = diff2[1] - d[1]*v[1][1];
   if(((*length2)=sqrt(diff2[0]*diff2[0]+diff2[1]*diff2[1])) > 0.1)
   {
      status = FALSE;
   }
   if(n > 50)
   {
      status = FALSE;
   }
   return(status);
}

/************************************************************
   Compute velocities
************************************************************/
void compute_vels(Ix,Iy,It,full_vels,norm_vels1,norm_vels2,E,pic_x,pic_y,n,tau_D,flag)
float Ix[PIC_X][PIC_Y],Iy[PIC_X][PIC_Y],It[PIC_X][PIC_Y],E[PIC_X][PIC_Y];
float full_vels[2][PIC_X][PIC_Y],norm_vels1[2][PIC_X][PIC_Y],tau_D;
float norm_vels2[2][PIC_X][PIC_Y];
int n,flag,pic_x,pic_y;
{
   float mag,M[2][2],MI[2][2],B[2],denominator;
   float eigenvalues[2],eigenvectors[2][2],length1,length2;
   float angle,temp1,v1,v2;
   float sigma1,sigma2,sigmap,temp,diff1[N],diff2[N];
   int i,j,k,l,full_count,norm_count1,norm_count2,no_count,nrot,eigen_count,no_swaps;
   float eigenvalues2[2],eigenvectors2[2][2],weight[5][5],sum,coeff[5],v[2];
   int mag_zero;
   int kk;


   printf("Eigenvalue Threshold: %f\n",tau_D);
   printf("Threshold on Raw Normal Velocities: %f\n",raw_mag);
   MAX_raw_mag = -HUGE;
   mag_zero = 0;
   fflush(stdout);
   /* Parameter values as specified in Simoncelli, Adelson and Heeger, page 313  */
   sigma1 = 0.08;
   sigma2 = 1.0;
   sigmap = 2.0;

   /* Compute weights */
   sum = 0.0;
   coeff[0] = coeff[4] = 0.0625;
   coeff[1] = coeff[3] = 0.25;
   coeff[2] = 0.375;
   for(i=0;i<5;i++)
   {
      for(j=0;j<5;j++)
      {
         weight[i][j] = coeff[i]*coeff[j];
         sum += weight[i][j];
      }
   }
   full_count = 0;
   norm_count1 = norm_count2 = 0;
   no_count = 0;
   eigen_count = 0;
   no_swaps = 0;
   fflush(stdout);
   for(i=0;i<PIC_X;i++)
      for(j=0;j<PIC_Y;j++)
      {
         full_vels[0][i][j] = full_vels[1][i][j] = NO_VALUE;
         norm_vels1[0][i][j] = norm_vels1[1][i][j] = NO_VALUE;
         norm_vels2[0][i][j] = norm_vels2[1][i][j] = NO_VALUE;
         /*if(DEBUG) for(kk=0;kk<8;kk++) data[kk][i][j] = 0.0;*/
         Imag[i][j] = E[i][j] = 0.0;
      }

   for(i=n;i<pic_x-n;i++)
      for(j=n;j<pic_y-n;j++)
      {
         M[0][0] = M[1][1] = M[0][1] = M[1][0] = 0.0;
         B[0] = B[1] = 0.0;
         mag = 1.0;
         /* Compute on 5*5 neighbourhood */
         for(k=(-2);k<2;k++)
            for(l=(-2);l<2;l++)
            {
               if(flag==TRUE) mag = sigma1*(Ix[i+k][j+l]*Ix[i+k][j+l]+
                   Iy[i+k][j+l]*Iy[i+k][j+l])+sigma2;
               M[0][0] = M[0][0] + weight[k+2][l+2]*(Ix[i+k][j+l]*Ix[i+k][j+l])/mag;
               M[1][1] = M[1][1] + weight[k+2][l+2]*(Iy[i+k][j+l]*Iy[i+k][j+l])/mag;
               M[0][1] = M[0][1] + weight[k+2][l+2]*(Ix[i+k][j+l]*Iy[i+k][j+l])/mag;
               B[0] = B[0] + weight[k+2][l+2]*(Ix[i+k][j+l]*It[i+k][j+l])/mag;
               B[1] = B[1] + weight[k+2][l+2]*(Iy[i+k][j+l]*It[i+k][j+l])/mag;
            }
         if(DEBUG)
         {
            data[0][i][j] = M[0][0];
            data[1][i][j] = M[0][1];
            data[2][i][j] = M[1][1];
            data[3][i][j] = B[0];
            data[4][i][j] = B[1];
            data[5][i][j] = Ix[i][j];
            data[6][i][j] = Iy[i][j];
            data[7][i][j] = It[i][j];
         }
         M[1][0] = M[0][1]; /* The M array is symmetric */
         if(flag==TRUE)
         {
            M[0][0] = M[0][0] + 1.0/sigmap;
            M[1][1] = M[1][1] + 1.0/sigmap;
         }
         /* Invert 2*2 matrix */
         denominator = M[0][0]*M[1][1]-M[1][0]*M[0][1]; /* The determinant of M */
         MI[0][0] = M[1][1]/denominator;
         MI[0][1] = -M[0][1]/denominator;
         MI[1][0] = -M[1][0]/denominator;
         MI[1][1] = M[0][0]/denominator;

         jacobi(M,2,eigenvalues,eigenvectors,&nrot);
         if(check_eigen_calc(M,eigenvalues,eigenvectors,nrot,diff1,diff2,
             &length1,&length2,&angle)==FALSE)
         {
            if(FALSE)
            {
               printf("\n********************************************\n");
               printf("Fatal error: eigenvalue/eigenvector error\n");
               printf("i=%d j=%d\n",i,j);
               printf("eigenvalues: %f %f\n",eigenvalues[0],eigenvalues[1]);
               printf("eigenvector1: %f %f\n",eigenvectors[0][0],eigenvectors[1][0]);
               printf("eigenvector2: %f %f\n",eigenvectors[0][1],eigenvectors[1][1]);
               printf("\n      M:                        MI\n");
               printf("%12.6f %12.6f %12.6f %12.6f\n",M[0][0],M[0][1],MI[0][0],MI[0][1]);
               printf("%12.6f %12.6f %12.6f %12.6f\n",M[1][0],M[1][1],MI[1][0],MI[1][1]);
               printf("B: %f %f\n",B[0],B[1]);
               printf("Determinant of M: %f\n",denominator);
               printf("nrot: %d\n",nrot);
               printf("Angle between two eigenvectors: %f degrees\n",angle);
               printf("Difference length for eigenvector1\n");
               printf("Difference: %f %f Length: %f\n",diff1[0],diff1[1],length1);
               printf("Difference length for eigenvector2\n");
               printf("Difference: %f %f Length: %f\n",diff2[0],diff2[1],length2);
               /* Check eigenvalues/eigenvectors for 2*2 matrix using
       closed form method as in Anandan's thesis */
               eigenvalues2[0] = 0.5*((M[0][0]+M[1][1]) - 
                   sqrt((M[0][0]-M[1][1])*(M[0][0]-M[1][1])+4.0*M[0][1]*M[1][0]));
               eigenvalues2[1] = 0.5*((M[0][0]+M[1][1]) +
                   sqrt((M[0][0]-M[1][1])*(M[0][0]-M[1][1])+4.0*M[0][1]*M[1][0]));
               angle = atan2(eigenvalues2[0]-M[0][0],M[0][1]);
               eigenvectors2[0][0] = -cos(angle);
               eigenvectors2[1][0] = -sin(angle);
               eigenvectors2[0][1] = -sin(angle);
               eigenvectors2[1][1] =  cos(angle);
               printf("\nUsing Anandan's calculation:\n");
               printf("Angle of rotation: %f degrees\n",angle*180/M_PI);
               printf("eigenvalues: %f %f\n",eigenvalues2[0],eigenvalues2[1]);
               printf("eigenvector1: %f %f\n",eigenvectors2[0][0],eigenvectors2[1][0]);
               printf("eigenvector2: %f %f\n",eigenvectors2[0][1],eigenvectors2[1][1]);
               check_eigen_calc(M,eigenvalues2,eigenvectors2,nrot,diff1,diff2,
                   &length1,&length2,&angle);
               printf("Angle between two eigenvectors: %f degrees\n",angle);
               printf("Difference length for eigenvector1\n");
               printf("Difference: %f %f Length: %f\n",diff1[0],diff1[1],length1);
               printf("Difference length for eigenvector2\n");
               printf("Difference: %f %f Length: %f\n",diff2[0],diff2[1],length2);
               printf("********************************************\n\n");
               fflush(stdout);
               if(FALSE) exit(1);
            }
            eigen_count++;
         }
         else 
         {
            /* Sort the eigenvalues and the corresponding eigenvectors */
            /* Most likely, already ordered           */
            if(eigenvalues[0] < eigenvalues[1]) /* Largest eigenvalue first */
            {
               /* swap eigenvalues */
               temp = eigenvalues[0];
               eigenvalues[0] = eigenvalues[1];
               eigenvalues[1] = temp;
               /* swap eigenvector components */
               temp = eigenvectors[0][0];
               eigenvectors[0][0] = eigenvectors[0][1];
               eigenvectors[0][1] = temp;
               temp = eigenvectors[1][0];
               eigenvectors[1][0] = eigenvectors[1][1];
               eigenvectors[1][1] = temp;
               no_swaps++;
            }

            /* Full velocity if spread of M is small */
            if(eigenvalues[0] >= tau_D && eigenvalues[1] >= tau_D)
            {
               if(denominator > 0.0)
               {
                  full_vels[1][i][j] = -(v1= -(MI[0][0]*B[0]+MI[0][1]*B[1]));
                  full_vels[0][i][j] =  (v2= -(MI[1][0]*B[0]+MI[1][1]*B[1]));
                  /* Compute the residual */
                  for(k=(-2);k<=2;k++)
                     for(l=(-2);l<=2;l++)
                     {
                        temp1 = weight[k+2][l+2]*
                            (Ix[i+k][j+l]*v1+
                            Iy[i+k][j+l]*v2+It[i+k][j+l]);
                        /* temp2 = weight[k+2][l+2]*It[i+k][j+l]; */
                        E[i][j] += (temp1*temp1);
                     }
                  full_count++;
               }
               else {
                  full_vels[0][i][j] = full_vels[1][i][j] = NO_VALUE;
               }
            }
            /* Normal velocity if spread of M is small in one direction only */
            else if(eigenvalues[0] > tau_D && fabs(denominator) > 0.00000001)
            /* Normal velocity if spread of MI is small in one direction only */
            {
               /* Project v onto that direction */
               v[0] = -(MI[0][0]*B[0]+MI[0][1]*B[1]);
               v[1] = -(MI[1][0]*B[0]+MI[1][1]*B[1]);
               norm_vels1[1][i][j] = -(v[0]*eigenvectors[0][0]+
                   v[1]*eigenvectors[1][0])*eigenvectors[0][0];
               norm_vels1[0][i][j] = (v[0]*eigenvectors[0][0]+
                   v[1]*eigenvectors[1][0])*eigenvectors[1][0];
               norm_count1++;
            }
            else 
            {
               no_count++;
            }
         }

         /* Compute type 2 normal velocity */
         mag = (Ix[i][j]*Ix[i][j] + Iy[i][j]*Iy[i][j]);
         Imag[i][j] = sqrt(mag);
         if(Imag[i][j] > MAX_raw_mag) MAX_raw_mag = Imag[i][j];
         if(Imag[i][j] > raw_mag)
         {
            norm_vels2[1][i][j] =  It[i][j]*Ix[i][j]/mag;
            norm_vels2[0][i][j] = -It[i][j]*Iy[i][j]/mag;
            norm_count2++;
         }
         else mag_zero++;
      }

   printf("%d full velocities computed\n",full_count);
   printf("%d least squares normal velocities computed\n",norm_count1);
   printf("%d raw normal velocities computed\n",norm_count2);
   printf("%d locations where velocity information thresholded\n",no_count);
   printf("%d locations where eigenvalue/eigenvector calculation failed\n",eigen_count);
   printf("%d locations with spatial gradient is zero\n",mag_zero);
   fflush(stdout);
   if(DEBUG) for(kk=0;kk<8;kk++) write(fd_temp[kk],&data[kk][0][0],150*150*4);
}


/************************************************************
   Output full velocities using old Burkitt format
************************************************************/
void output_velocities(fdf,s,full_velocities,pic_x,pic_y,n)
float full_velocities[2][PIC_X][PIC_Y];
char s[32];
int fdf,n;
{
   float x,y;
   int i,j,bytes,no_novals,no_vals,NORMAL;

   bytes = 0;
   NORMAL = FALSE;
   if(strcmp(s,"Normal")==0) NORMAL = TRUE;
   if(fdf==NULL)
   {
      printf("\nFatal error: full velocity file not opened\n");
      exit(1);
   }
   /* original size */
   y = pic_x;
   x = pic_y;
   bytes += write(fdf,&x,sizeof(int));
   bytes += write(fdf,&y,sizeof(int));

   /* size of result data */
   y = pic_x-2*n;
   x = pic_y-2*n;
   bytes += write(fdf,&x,sizeof(int));
   bytes += write(fdf,&y,sizeof(int));

   /* offset to start of data */
   y = n;
   x = n;
   bytes += write(fdf,&x,sizeof(int));
   bytes += write(fdf,&y,sizeof(int));

   no_novals = no_vals = 0;
   /* Prepare velocities for output, i.e. rotate by 90 degrees */
   for(i=n;i<pic_x-n;i++)
      for(j=n;j<pic_y-n;j++)
      {
         if(full_velocities[0][i][j] != NO_VALUE && 
             full_velocities[1][i][j] != NO_VALUE)
         {
            no_vals++;
         }
         else
         {
            no_novals++;
         }
      }
   for(i=n;i<pic_x-n;i++)
      for(j=n;j<pic_y-n;j++)
   {
      bytes += write(fdf,&full_velocities[0][i][j],sizeof(float));
      bytes += write(fdf,&full_velocities[1][i][j],sizeof(float));
   }
   close(fdf);
   printf("\n%s velocities output from output_velocities: %d bytes\n",s,bytes);
   printf("Number of positions with velocity: %d\n",no_vals);
   printf("Number of positions without velocity: %d\n",no_novals);
   printf("Percentage of %s velocities: %f\n",s,
       no_vals/(1.0*(no_vals+no_novals))*100.0);
   fflush(stdout);
}

/*****************************************************************/
/* Read the correct velocity data                                */
/*****************************************************************/
void read_correct_vels (fname, correct_vels)
char fname[];
float correct_vels[2][PIC_X][PIC_Y];
{
   int i, fd, no_bytes=0;

   if ((fd = open(fname,O_RDONLY))==(-1))
   {
      printf("Fatal error in opening file %s\n", fname);
      printf("fd_correct: %d\n",fd);
      exit(1);
   }
   no_bytes += read (fd, &actual_y, 4);
   no_bytes += read (fd, &actual_x, 4);
   no_bytes += read (fd, &size_y, 4);
   no_bytes += read (fd, &size_x, 4);
   no_bytes += read (fd, &offset_y, 4);
   no_bytes += read (fd, &offset_x, 4);
   if (offset_x != 0.0 || offset_y != 0.0 ||
       actual_x != size_x || actual_y != size_y)
   {
      printf ("Fatal error: something wrong with correct velocity data\n");
      printf ("Actual y: %f Actual x: %f\n",actual_y,actual_x);
      printf ("Size y: %f Size x: %f\n",size_y,size_x);
      printf ("Offset y: %f Offset x: %f\n",offset_y,offset_x);
      exit (1);
   }
   int_size_y = size_y;
   int_size_x = size_x;
   for (i=0;i<int_size_x;i++)
      no_bytes += read (fd, &correct_vels[0][i][0], int_size_y*8);
   printf ("\nFile %s opened and read\n",fname);
   printf ("Size of correct velocity data: %d %d\n",int_size_y,int_size_x);
   printf ("%d bytes read\n",no_bytes);
   fflush (stdout);
   close (fd);
}


/***************************************************************/
/*  Compute error statistics                                   */
/***************************************************************/

/************************************************************
   norm of a vector v of length n.
************************************************************/
float norm(v,n)
float v[];
int n;
{
   int i;
   float sum = 0.0;

   for (i=0;i<n; i++)
      sum += (v[i]*v[i]);
   sum = sqrt(sum);
   return sum;
}

/************************************************************
 Full Image Velocity Angle Error
************************************************************/
float PsiER(ve,va)
float ve[2],va[2];
{
   float nva, nve;
   float v,r;
   float VE[3],VA[3];

   VE[0] = ve[0];
   VE[1] = ve[1];
   VE[2] = 1.0;

   VA[0] = va[0];
   VA[1] = va[1];
   VA[2] = 1.0;

   nva = norm(VA,3);
   nve = norm(VE,3);
   v = (VE[0]*VA[0]+VE[1]*VA[1]+1.0)/(nva*nve);

   /**  sometimes roundoff error causes problems **/
   if(v>1.0 && v < 1.0001) v = 1.0;

   r = acos(v)*180.0/PI;
   if (!(r>=0.0 && r< 180.0))
   {
      printf("ERROR in PSIER()...\n r=%8.4f v=%8.4f nva=%8.4f nve= %8.4f\n",
          r,v,nva,nve);
      printf("va=(%f,%f) ve=(%f,%f)\n",va[0],va[1],ve[0],ve[1]);
   }
   return r;
}


/************************************************************
 Normal Image Velocity Angle Error
************************************************************/
float PsiEN(ve,va)
float ve[2],va[2];
{
   float nva,nve;
   float v1,v2;
   float n[2];

   nva = norm(va,2), nve = norm(ve,2);
   if(nve > 0.00000001)
   {
      n[0] = ve[0]/nve;
      n[1] = ve[1]/nve;
      v1 = (va[0]*n[0] + va[1]*n[1]-nve) ;
      v2 = v1/(sqrt((1.0+nva*nva))*sqrt((1.0+nve*nve)));
      v1 =  asin(v2)*180.0/PI;
      if(!(v1>=-90.0 && v1<=90.0))
      {
         printf("ERROR in PSIEN()  v1: %f ve: %f\n",v1,v2);
         printf("nve: %f nva: %f\n",nve,nva);
         printf("n: %f %f\n",n[0],n[1]);
         printf(" ve: %f %f va: %f %f\n",ve[0],ve[1],va[0],va[1]);
         fflush(stdout);
      }
   }
   else v1 = NO_VALUE;
   return fabs(v1);
}


/************************************************************
    Compute error statistics                                
************************************************************/
void calc_statistics(correct_fname,norm_vels1,norm_vels2,int_size_x,int_size_y,full_vels,E,
pic_x,pic_y,n,ave_error,st_dev,density,residual,min_angle,max_angle,
norm_ave_error1,norm_st_dev1,norm_density1,norm_min_angle1,norm_max_angle1,
norm_ave_error2,norm_st_dev2,norm_density2,norm_min_angle2,norm_max_angle2)
char correct_fname[32];
float full_vels[2][PIC_X][PIC_Y],*ave_error,*density,*st_dev,*residual;
float E[PIC_X][PIC_Y],*min_angle,*max_angle;
float norm_vels1[2][PIC_X][PIC_Y],norm_vels2[2][PIC_X][PIC_Y];
float *norm_min_angle1,*norm_max_angle1,*norm_ave_error1,*norm_st_dev1,*norm_density1;
float *norm_min_angle2,*norm_max_angle2,*norm_ave_error2,*norm_st_dev2,*norm_density2;
int n,pic_x,pic_y,int_size_x,int_size_y;
{
   int full_count,norm_count1,no_full_count,no_norm_count1,i,j,total_count;
   int norm_count2,no_norm_count2;
   float normal_sumX2_1,normal_sumX2_2,sumX2,temp,uva[2],uve[2],bin_err2;
   float histogram_error[100],histogram_error2[100],bin_sum,bin_ave,bin_st_dev;
   float bin_density;
   int histogram_count[100],int_mag,bin_ct;
   float correct_vels[2][PIC_X][PIC_Y];

   full_count = norm_count1 = norm_count2 = no_full_count = no_norm_count1 = no_norm_count2 = total_count = 0;
   normal_sumX2_1 = normal_sumX2_2 = sumX2 = 0.0;
   (*norm_min_angle1) = (*min_angle) = HUGE;
   (*norm_max_angle1) = (*max_angle) = 0.0;
   (*norm_min_angle2) = HUGE;
   (*norm_max_angle2) = 0.0;
   (*ave_error) = (*st_dev) = (*density) = (*residual) = 0.0;
   (*norm_ave_error1) = (*norm_st_dev1) = (*norm_density1) = 0.0;
   (*norm_ave_error2) = (*norm_st_dev2) = (*norm_density2) = 0.0;

   for(i=0;i<100;i++) {
      histogram_error[i]=histogram_error2[i]=0.0;
      histogram_count[i]=0;
   }

   read_correct_vels (correct_fname, correct_vels);

   for(i=n;i<pic_x-n;i++)
   {
      for(j=n;j<pic_y-n;j++)
      {
         if(full_vels[0][i][j] != NO_VALUE && full_vels[1][i][j] != NO_VALUE)
         {
            full_count++;
            uve[0] = full_vels[0][i][j];
            uve[1] = full_vels[1][i][j];
            uva[0] = correct_vels[0][i][j];
            uva[1] = correct_vels[1][i][j];
            temp = PsiER(uve,uva);
            (*ave_error) += temp;
            if(E[i][j] == NO_VALUE)
            {
               printf("Fatal error: E has no value\n");
               exit(1);
            }
            (*residual) += E[i][j];
            sumX2 += temp*temp;
            if(temp < (*min_angle)) (*min_angle) = temp;
            if(temp > (*max_angle)) (*max_angle) = temp;
         }
         else no_full_count++;

         if(norm_vels1[0][i][j] != NO_VALUE && norm_vels1[1][i][j] != NO_VALUE)
         {
            norm_count1++;
            uve[0] = norm_vels1[0][i][j];
            uve[1] = norm_vels1[1][i][j];
            uva[0] = correct_vels[0][i][j];
            uva[1] = correct_vels[1][i][j];
            temp = PsiEN(uve,uva);
            (*norm_ave_error1) += temp;
            normal_sumX2_1 += temp*temp;
            if(temp < (*norm_min_angle1)) (*norm_min_angle1) = temp;
            if(temp > (*norm_max_angle1)) (*norm_max_angle1) = temp;
         }
         else no_norm_count1++;

         if(norm_vels2[0][i][j] != NO_VALUE && norm_vels2[1][i][j] != NO_VALUE)
         {
            norm_count2++;
            uve[0] = norm_vels2[0][i][j];
            uve[1] = norm_vels2[1][i][j];
            uva[0] = correct_vels[0][i][j];
            uva[1] = correct_vels[1][i][j];
            temp = PsiEN(uve,uva);
            (*norm_ave_error2) += temp;
            normal_sumX2_2 += temp*temp;
            if(temp < (*norm_min_angle2)) (*norm_min_angle2) = temp;
            if(temp > (*norm_max_angle2)) (*norm_max_angle2) = temp;
            int_mag = (Imag[i][j]);
            if(int_mag >= MAX_I) int_mag = MAX_I;
            histogram_error[int_mag] += temp;
            histogram_error2[int_mag] += temp*temp;
            histogram_count[int_mag]++;
         }
         else no_norm_count2++;

         total_count++;
      }
   }
   (*density) = (full_count*100.0)/(total_count);
   (*norm_density1) = (norm_count1*100.0)/(total_count);
   (*norm_density2) = (norm_count2*100.0)/(total_count);

   if(full_count != 0) (*ave_error) = (*ave_error)/full_count;
   else (*ave_error) = 0.0;
   if(norm_count1 != 0) (*norm_ave_error1) = (*norm_ave_error1)/norm_count1;
   else (*norm_ave_error1) = 0.0;
   if(norm_count2 != 0) (*norm_ave_error2) = (*norm_ave_error2)/norm_count2;
   else (*norm_ave_error2) = 0.0;

   if(full_count > 1)
   {
      temp = fabs((sumX2 - full_count*(*ave_error)*(*ave_error))/(full_count-1));
      (*st_dev) = sqrt(temp);
   }
   else (*st_dev) = 0.0;
   if(norm_count1 > 1)
      (*norm_st_dev1) = sqrt((normal_sumX2_1 - norm_count1*(*norm_ave_error1)*(*norm_ave_error1))/(norm_count1-
          1));
   else (*norm_st_dev1) = 0.0;
   if(norm_count2 > 1)
      (*norm_st_dev2) = sqrt((normal_sumX2_2 - norm_count2*(*norm_ave_error2)*(*norm_ave_error2))/(norm_count2-
          1));
   else (*norm_st_dev2) = 0.0;

   if(full_count != 0) (*residual) = (*residual)/full_count;
   if((*ave_error) == 0.0) {
      (*min_angle) = (*max_angle) = 0.0;
   }
   if((*norm_ave_error1) == 0.0) {
      (*norm_min_angle1) = (*norm_max_angle1) = 0.0;
   }
   if((*norm_ave_error2) == 0.0) {
      (*norm_min_angle2) = (*norm_max_angle2) = 0.0;
   }

   printf("\nIn calc_statistics\n");
   printf("%d full velocities\n",full_count);
   printf("%d positons without full velocity\n",no_full_count);
   printf("%d positions in total\n",total_count);
   printf("%d least squares normal velocities\n",norm_count1);
   printf("%d positons without least squares normal velocity\n",no_norm_count1);
   printf("%d positions with full or least squares normal velocity\n",
       total_count-full_count-norm_count1);
   printf("%d raw normal velocities\n",norm_count2);
   printf("%d positons without raw normal velocity\n",no_norm_count2);

   fp =fopen("lucas.plot.data","w");
   fprintf(fp,"%d\n",MAX_I);
   if(raw_statistics)
   {
      printf("\n                             Raw Normal Velocity Analysis\n");
      printf("\n                                  Histogram Error\n\n");
      for(i=0;i<MAX_I;i++)
      {
         if(histogram_count[i] != 0.0) bin_ave = histogram_error[i]/(histogram_count[i]*1.0);
         else bin_ave = 0.0;
         if(histogram_count[i] > 1) temp = ((histogram_error2[i]-histogram_count[i]*bin_ave*bin_ave)/(histogram_count[i]-
             1.0));
         else temp = 0.0;
         bin_density = histogram_count[i]*100.0/total_count;
         bin_st_dev = sqrt(fabs(temp));
         printf("Bin:%5.2f Average Error:%8.5f St. Dev:%8.5f Count:%5d Density:%6.3f\n",i+0.5,bin_ave,bin_st_dev,
             histogram_count[i],bin_density);
         fprintf(fp,"%f %f %f %f\n",i+0.5,bin_ave,bin_st_dev,bin_density);
      }
      fprintf(fp,"\n\n\n");
      fprintf(fp,"%d\n",MAX_I);
      printf("\n                                 Cumulative Error\n\n");
      bin_sum = 0.0;
      bin_ct=0;
      bin_err2 = 0.0;
      for(i=0;i<MAX_I;i++)
      {
         bin_sum += histogram_error[i];
         bin_err2 += histogram_error2[i];
         bin_ct += histogram_count[i];
      }
      for(i=0;i<MAX_I;i++)
      {
         if(bin_ct != 0) bin_ave = bin_sum/bin_ct;
         else bin_ave = 0.0;
         bin_density = bin_ct*100.0/total_count;
         if(bin_ct > 1)  temp = ((bin_err2-bin_ct*bin_ave*bin_ave)/(bin_ct-1));
         else temp = 0.0;
         bin_st_dev = sqrt(fabs(temp));
         printf("Bin:>%5.2f Average Error:%8.5f St. Dev:%8.5f Count:%5d Density:%6.3f\n",i*1.0,bin_ave,bin_st_dev,
             bin_ct,bin_density);
         fprintf(fp,"%f %f %f %f\n",i*1.0,bin_ave,bin_st_dev,bin_density);
         bin_sum -= histogram_error[i];
         bin_err2 -= histogram_error2[i];
         bin_ct -= histogram_count[i];
      }
      fprintf(fp,"\n");
      fclose(fp);
   }
   fflush(stdout);
}

/*********************************************************************/
/*   Main program                  */
/*********************************************************************/

int main (argc, argv)
int argc;
char **argv;
{
   unsigned char header[HEAD];
   char path1[32], path2[32], path3[32];
   char correct_filename[32];
   char full_name[64], norm_name[64], raw_name[64];
   int fdf, fdn, fdr, fd_correct;
   int size, offset, start, middle, end, no_bytes;
   int flag, i;
   float sigma, tau_D;
   float ave_error,st_dev,residual,density,min_angle,max_angle;
   float norm_ave_error1,norm_st_dev1,norm_density1,
   norm_min_angle1,norm_max_angle1;
   float norm_ave_error2,norm_st_dev2,norm_density2,
   norm_min_angle2,norm_max_angle2;

   if(argc < 6 || argc > 17)
   {
      usage(argv[0]);
      exit(1);
   }
   strcpy(path1,".");
   strcpy(path2,".");
   strcpy(path3,".");
   sscanf(argv[2],"%f",&sigma);
   sscanf(argv[3],"%d",&middle);
   sscanf(argv[4],"%f",&tau_D);
   strcpy(path1,argv[5]);
   strcpy(path3,argv[6]);
   /* Show command line arguments */
   printf("\n--------------------------------\nCommand line: ");
   for(i=0;i<argc;i++)
      printf("%s ",argv[i]);
   printf("\n%d arguments\n",argc-1);
   printf("sigma=%f\n",sigma);
   printf("Central image: %d\n",middle);
   printf("Tau threshold: %f\n",tau_D);
   printf("Input directory: %s\n",path1);
   printf("Output directory: %s\n",path3);
   fflush(stdout);

   binary = FALSE;
   output_smooth = FALSE;
   raw_statistics = FALSE;
   raw_mag = 0.0;
   flag = FALSE;
   i = 7;
   strcpy(correct_filename,"unknown");
   while(i<argc)
   {
      if(strcmp("-M",argv[i])==0)
      {
         flag = TRUE;
         i++;
      }
      else if(strcmp("-L",argv[i])==0) i++;
      else if(strcmp("-C",argv[i])==0)
      {
         strcpy(correct_filename,argv[i+1]);
         i+=2;
      }
      else if(strcmp("-T",argv[i])==0)
      {
         sscanf(argv[i+1],"%f",&raw_mag);
         if(raw_mag==0.0) raw_statistics = TRUE;
         i += 2;
      }
      else if(strcmp("-B",argv[i])==0)
      {
         sscanf(argv[i+1],"%d",&pic_y);
         sscanf(argv[i+2],"%d",&pic_x);
         binary = TRUE;
         i += 3;
      }
      else if(strcmp("-S",argv[i])==0)
      {
         strcpy(path2,argv[i+1]);
         printf("Smoothed data directory: %s\n",path2);
         output_smooth = TRUE;
         i += 2;
      }
      else
      {
         printf("Invalid option %s - program terminates\n",argv[i]);
         exit(1);
      }
   }
   printf("Correct velocity file: %s\n",correct_filename);
   fflush(stdout);

   size = 6*sigma+1;
   if(size%2==0) size = size+1;
   offset = size/2+2; /* Add 2 as neighbourhood size offset */
   start = middle-offset;
   end = middle+offset;
   printf("Size: %d Offset: %d Start: %d End: %d\n",size,offset,start,end);
   printf("%d images required\n",end-start+1);
   if(end < start)
   {
      printf("Specify images in ascending order\n");
      exit(1);
   }

   Smooth3D_Gaussian (path1, argv[1], floatpic, sigma, pic_t, pic_x, pic_y,
                      start, middle, end, header);
   if(output_smooth && sigma!= 0.0)
      writefiles (path2, argv[1], floatpic, sigma,
          pic_t, pic_x, pic_y,
          middle-2, middle+2,
          header);
   else if(sigma != 0.0)
      printf("\nSmoothed images not output\n");

   if(flag==FALSE)
   {
      sprintf(full_name,"%s/lucas.%sF-%4.4f-%3.1f",path3,argv[1],tau_D,sigma);
      sprintf(norm_name,"%s/lucas.%sN-%4.4f-%3.1f",path3,argv[1],tau_D,sigma);
      sprintf(raw_name,"%s/lucas.%sR-%4.4f-%3.1f",path3,argv[1],tau_D,sigma);
   }
   else
   {
      sprintf(full_name,"%s/Mlucas.%sF-%4.4f-%3.1f",path3,argv[1],tau_D,sigma);
      sprintf(norm_name,"%s/Mlucas.%sN-%4.4f-%3.1f",path3,argv[1],tau_D,sigma);
      sprintf(raw_name,"%s/Mlucas.%sR-%4.4f-%3.1f",path3,argv[1],tau_D,sigma);
   }
   printf("Output full velocities go to file %s\n",full_name);
   printf("Output least squares normal velocities go to file %s\n",norm_name);
   printf("Output raw normal velocities go to file %s\n",raw_name);
   printf("\n");
   fflush(stdout);

   compute_ders (floatpic, Ix, Iy, It, pic_t, pic_x, pic_y, offset);
   compute_vels (Ix, Iy, It, full_vels, &floatpic[0], &floatpic[2],
       &floatpic[4], pic_x, pic_y, 2*offset, tau_D, flag);
   if((fdf=creat(full_name,0644))!=NULL)
      output_velocities(fdf,"Full",full_vels,pic_x,pic_y,2*offset);
   else
      printf("Error in opening %s file\n\n",full_name);
   if((fdn=creat(norm_name,0644))!=NULL)
      output_velocities(fdn,"Normal",&floatpic[0],pic_x,pic_y,2*offset);
   else
      printf("Error in opening %s file\n\n",norm_name);
   if((fdr=creat(raw_name,0644))!=NULL)
      output_velocities(fdn,"Normal",&floatpic[2],pic_x,pic_y,2*offset);
   else printf("Error in opening %s file\n\n",raw_name);
   printf("\n");
   fflush(stdout);

   if(strcmp(correct_filename,"unknown")!=0)
   {
      calc_statistics(correct_filename, &floatpic[0], &floatpic[2],
          int_size_x, int_size_y, full_vels, &floatpic[4],
          pic_x, pic_y, 2*offset,
          &ave_error, &st_dev, &density, &residual,
          &min_angle, &max_angle,
          &norm_ave_error1, &norm_st_dev1,&norm_density1,
          &norm_min_angle1, &norm_max_angle1,
          &norm_ave_error2, &norm_st_dev2, &norm_density2,
          &norm_min_angle2,&norm_max_angle2);
      printf("\n\nComputed Statistics\n");
      printf("Lambda_2: %f  Error: %f St Dev: %f\n",tau_D,ave_error,st_dev);
      printf("Lambda_2: %f  Density: %f\n",tau_D,density);
      printf("Lambda_2: %f  Residual: %f\n",tau_D,residual);
      printf("Minimum angle error: %f Maximum angle error: %f\n",
          min_angle,max_angle);
      printf("\nLeast Squares Normal Velocity Results\n");
      printf("Normal Error: %f Normal Standard Deviation: %f\n",
          norm_ave_error1,norm_st_dev1);
      printf("Normal Density: %f\nMinimum Normal Angle: %f Maximum Normal Angle: %f\n",
          norm_density1,norm_min_angle1,norm_max_angle1);
      printf("\nRaw Normal Velocity Results\n");
      printf("Normal Error: %f Normal Standard Deviation: %f\n",
          norm_ave_error2,norm_st_dev2);
      printf("Normal Density: %f\nMinimum Normal Angle: %f Maximum Normal Angle: %f\n",
          norm_density2,norm_min_angle2,norm_max_angle2);
   }
   fflush(stdout);
   exit (0);
}
