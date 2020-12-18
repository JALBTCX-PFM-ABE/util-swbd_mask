
/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.
*********************************************************************************************/

/****************************************  IMPORTANT NOTE  **********************************

    Comments in this file that start with / * ! or / / ! are being used by Doxygen to
    document the software.  Dashes in these comment blocks are used to create bullet lists.
    The lack of blank lines after a block of dash preceeded comments means that the next
    block of dash preceeded comments is a new, indented bullet list.  I've tried to keep the
    Doxygen formatting to a minimum but there are some other items (like <br> and <pre>)
    that need to be left alone.  If you see a comment that starts with / * ! or / / ! and
    there is something that looks a bit weird it is probably due to some arcane Doxygen
    syntax.  Be very careful modifying blocks of Doxygen comments.

*****************************************  IMPORTANT NOTE  **********************************/


/*******************************************************************************************/
/*!

  - Module Name:        swbd_mask

  - Programmer(s):      Jan C. Depner (PFM Software)

  - Date Written:       October 2013

  - Purpose:            Reads the uncompressed SRTM Water Body Data (SWBD) shape files and
                        creates a world (or as much as is covered) land mask.

  - Arguments:          argv[1]         -   resolution (in seconds - 1, 3, 10, 30, or 60)

  - Caveats:            You must have all of the uncompressed SWBD files in a single
                        directory in order to run this.  You must also have the dummy
                        *.wtr and *.lnd files for the cells that don't have associated
                        shape files.  Make sure that the ABE_DATA environment variable points
                        to the directory that holds the SWBD directory and the land_mask
                        directory and that you have write access to the land_mask directory
                        since that is where the output file will be built.  The output file
                        will be named swbd_mask_XX_second.clm, where XX is 01, 03, 10, 30,
                        or 60.


  - Description of the compressed land mask (.clm) file format (look Ma, no endians!)


    <pre>

        Header - 16384 bytes, ASCII

        [HEADER SIZE] = 16384
        [CREATION DATE] = 
        [VERSION] = 
        [ZLIB VERSION] =
        [RESOLUTION] = 1, 3, 10, 30, or 60
        [END OF HEADER]


        One-degree map - 64800 * 7 bytes, binary, stored as unsigned characters.

            Single record (7 bytes) :

                32 bits - 0 = undefined, 1 = all land, 2 = all water, otherwise the address
                          of the compressed block
                24 bits - 0 or the size of the compressed block (*CBS) 

            Records start at 90S,180W and proceed west to east then south to north (that is,
            the second record is for 90S,179W and the 361st record is for 89S,180W).


        Data - 1's and 0's (woo hoo)

            *CBS bytes - data

            The data is stored as a series of single bits for water (0) and land (1).  Each
            bit represents a one, three, ten, thirty, or sixty second cell in the block.  The
            block is a one-degree square.  It will be 3600 X 3600, 1200 X 1200, 360 X 360,
            120 X 120, or 60 X 60 depending on the resolution.  It is ordered west to east
            starting in the southwest corner and moving northward.  The compression is
            compliments of the ZLIB compression library which can be found at
            http://www.zlib.net/.  Many thanks to Jean-loup Gailly, Mark Adler, and all
            others associated with that effort.

    </pre>

*********************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <zlib.h>

#include "nvutility.hpp"
#include "nvutility.h"

#include "shapefil.h"
#include "version.h"

#include "maskThread.hpp"


#define SWBD_MASK_HEADER_SIZE        16384


void usage (char *string)
{
  fprintf (stderr, "Usage: %s RESOLUTION [NUM_THREADS]\n\n", string);
  fprintf (stderr, "Where\n");
  fprintf (stderr, "\tRESOLUTION = resolution of mask in seconds (1, 3, 10, 30, or 60)\n");
  fprintf (stderr, "\tNUM_THREADS = number of compute threads (4[default] or 16)\n\n");
  exit (-1);
}



int32_t main (int32_t argc, char **argv)
{
  int32_t           type, numShapes, resolution = 0, total_blocks = 0, num_threads = 4;
  double            minBounds[4], maxBounds[4], total_block_size = 0.0;
  char              dirname[512], shpname[512], lathem, lonhem, dataset[7] = {'a', 'e', 'f', 'i', 'n', 's', 'x'};
  FILE              *ofp;
  char              ofile[512];
  maskThread        mask_thread[16];


  printf ("\n\n%s\n\n", VERSION);


  if (argc < 2) usage (argv[0]);


  //  Check for ABE_DATA environment variable.

  if (getenv ("ABE_DATA") == NULL)
    {
      fprintf (stderr, "\n\nEnvironment variable ABE_DATA is not set\n\n");
      fflush (stderr);
      exit (-1);
    }

  sprintf (dirname, "%s", getenv ("ABE_DATA"));


  if (check_srtm_mask (3))
    {
      fprintf (stderr, "Can't find 3 second SRTM landmask\n\n");
      exit (-1);
    }


  sscanf (argv[1], "%d", &resolution);

  if (resolution != 1 && resolution != 3 && resolution != 10 && resolution != 30 && resolution != 60) usage (argv[0]);


  if (argc == 3) sscanf (argv[2], "%d", &num_threads);


  if (num_threads != 4 && num_threads != 16) usage (argv[1]);


  //  Open the output file.

  sprintf (ofile, "%s%1cland_mask%1cswbd_mask_%02d_second.clm", dirname, (char) SEPARATOR, (char) SEPARATOR, resolution);

  if ((ofp = fopen (ofile, "wb+")) == NULL)
    {
      perror (ofile);
      exit (-1);
    }


  //  Write the (minimalist) ASCII header.

  time_t t = time (&t);
  struct tm *cur_tm = gmtime (&t);

  fprintf (ofp, "[HEADER SIZE] = %d\n", SWBD_MASK_HEADER_SIZE);
  fprintf (ofp, "[VERSION] = %s\n", VERSION);
  fprintf (ofp, "[ZLIB VERSION] = %s\n", zlibVersion ());
  fprintf (ofp, "[CREATION DATE] = %s", asctime (cur_tm));
  fprintf (ofp, "[RESOLUTION] = %d\n", resolution);
  fprintf (ofp, "[END OF HEADER]\n");


  //  Zero out the remainder of the header.

  uint8_t zero = 0;
  int32_t j = ftell (ofp);
  for (int32_t i = j ; i < SWBD_MASK_HEADER_SIZE ; i++) fwrite (&zero, 1, 1, ofp);


  //  Set the default for all map addresses to 0 (undefined).

  uint8_t mapbuf[7];
  memset (mapbuf, 0, 7);

  for (int32_t i = -90 ; i < 90 ; i++)
    {
      for (int32_t j = -180 ; j < 180 ; j++)
        {
          fwrite (mapbuf, 7, 1, ofp);
        }
    }


  double **poly_x = NULL;
  double **poly_y = NULL;
  int32_t *poly_count = NULL;
  uint8_t *block = NULL;
  uint8_t *bit_block = NULL;
  SHPHandle shpHandle;
  SHPObject *shape = NULL;
  FILE *tfp = NULL;

  for (int32_t lat = -90 ; lat < 90 ; lat++)
    {
      if (lat < 0)
        {
          lathem = 's';
        }
      else
        {
          lathem = 'n';
        }

      int32_t lt = abs (lat);


      for (int32_t lon = -180 ; lon < 180 ; lon++)
        {
          //if (lat == 10 && lon == 124){
          uint32_t map_address = SWBD_MASK_HEADER_SIZE + ((lat + 90) * 360 + (lon +180)) * 7;


          if (lon < 0)
            {
              lonhem = 'w';
            }
          else
            {
              lonhem = 'e';
            }


          int32_t ln = abs (lon);


          //  Initialize variables

          int32_t num_poly = -1;
          uint8_t file_check = NVFalse;
          int32_t file_ext = 0;


          //  Check to make sure we have a valid file before we open the output file.

          for (int32_t ds = 0 ; ds < 6 ; ds++)
            {
              sprintf (shpname, "%s%1cSWBD%1c%1c%03d%1c%02d%1c.shp", dirname, (char) SEPARATOR, (char) SEPARATOR, lonhem, ln, lathem, lt, dataset[ds]);


              //  Make sure the file exists before we try to open it with the shape library.

              if ((tfp = fopen (shpname, "rb")) != NULL)
                {
                  file_check = NVTrue;
                  file_ext = ds;
                  fclose (tfp);
                  break;
                }
            }


          //  If we didn't find a file, check SRTM3 for all water or land, or, if we're outside of 57S to 60N, set to undefined.

          if (!file_check)
            {
              //  Undefined

              if (lat < -57 || lat > 59)
                {
                  memset (mapbuf, 0, 7);
                  fseek (ofp, map_address, SEEK_SET);
                  fwrite (mapbuf, 7, 1, ofp);
                }
              else
                {
                  double slat = lat + 0.5;
                  double slon = lon + 0.5;

                  int32_t lnd = read_srtm_mask_min_res (slat, slon, 3);


                  //  Water

                  if (lnd == 0)
                    {
                      memset (mapbuf, 0, 7);
                      bit_pack (mapbuf, 0, 32, 2);
                      fseek (ofp, map_address, SEEK_SET);
                      fwrite (mapbuf, 7, 1, ofp);
                    }


                  //  Land

                  else
                    {
                      memset (mapbuf, 0, 7);
                      bit_pack (mapbuf, 0, 32, 1);
                      fseek (ofp, map_address, SEEK_SET);
                      fwrite (mapbuf, 7, 1, ofp);
                    }
                }
            }
          else
            {
              //  If we have a shape file, define the name and process it.

              sprintf (shpname, "%s%1cSWBD%1c%1c%03d%1c%02d%1c.shp", dirname, (char) SEPARATOR, (char) SEPARATOR, lonhem, ln, lathem, lt, dataset[file_ext]);


              //  Open shape file

              shpHandle = SHPOpen (shpname, "rb");

              if (shpHandle == NULL)
                {
                  perror (shpname);
                  exit (-1);
                }


              fprintf (stderr,"Reading %s                        \n", shpname);
              fflush (stderr);


              //  Get shape file header info

              SHPGetInfo (shpHandle, &numShapes, &type, minBounds, maxBounds);


              //  Read all shapes

              for (int32_t i = 0 ; i < numShapes ; i++)
                {
                  shape = SHPReadObject (shpHandle, i);


                  //  Get all vertices

                  if (shape->nVertices >= 2)
                    {
                      for (int32_t j = 0, numParts = 1 ; j < shape->nVertices ; j++)
                        {
                          uint8_t start_segment = NVFalse;


                          //  Check for start of a new segment.

                          if (!j && shape->nParts > 0) start_segment = NVTrue;


                          //  Check for the start of a new segment inside a larger group of points (this would be a "Ring" point).

                          if (numParts < shape->nParts && shape->panPartStart[numParts] == j)
                            {
                              start_segment = NVTrue;
                              numParts++;
                            }


                          //  Start a new segment

                          if (start_segment)
                            {
                              //  Since num_poly starts at -1 this is perfectly cool.

                              num_poly++;


                              //  Allocate the count array.

                              poly_count = (int32_t *) realloc (poly_count, (num_poly + 1) * sizeof (int32_t));
                              if (poly_count == NULL)
                                {
                                  perror ("Allocating poly_count memory");
                                  exit (-1);
                                }


                              //  Set the count for the new arrays to zero.

                              poly_count[num_poly] = 0;


                              //  Allocate the polygon arrays.

                              poly_x = (double **) realloc (poly_x, (num_poly + 1) * sizeof (double *));
                              if (poly_x == NULL)
                                {
                                  perror ("Allocating poly_x memory");
                                  exit (-1);
                                }
                              poly_x[num_poly] = NULL;


                              poly_y = (double **) realloc (poly_y, (num_poly + 1) * sizeof (double *));
                              if (poly_y == NULL)
                                {
                                  perror ("Allocating poly_y memory");
                                  exit (-1);
                                }
                              poly_y[num_poly] = NULL;
                            }


                          //  Allocate memory for the new point.

                          poly_x[num_poly] = (double *) realloc (poly_x[num_poly], (poly_count[num_poly] + 1) * sizeof (double));
                          if (poly_x[num_poly] == NULL)
                            {
                              perror ("Allocating poly_x[num_poly] memory");
                              exit (-1);
                            }

                          poly_y[num_poly] = (double *) realloc (poly_y[num_poly], (poly_count[num_poly] + 1) * sizeof (double));
                          if (poly_y[num_poly] == NULL)
                            {
                              perror ("Allocating poly_y[num_poly] memory");
                              exit (-1);
                            }


                          //  Add point to current segment

                          poly_x[num_poly][poly_count[num_poly]] = shape->padfX[j];
                          poly_y[num_poly][poly_count[num_poly]] = shape->padfY[j];


                          //  Increment the point counter.

                          poly_count[num_poly]++;
                        }
                    }


                  //  Destroy the shape object.

                  SHPDestroyObject (shape);
                }


              //  Increment num_poly to account for the last polygon.

              num_poly++;


              //  Close the input file.

              SHPClose (shpHandle);


              //  Allocate the uint8_t block for the threads to put the land/water flags into.
              //  We have to use a block that is byte aligned so that the threads don't step
              //  on each other (as could happen if we tried to use bit_pack to set bits in
              //  a bit block).

              int32_t point_count = 3600 / resolution;

              block = (uint8_t *) calloc (point_count * point_count, sizeof (uint8_t));

              if (block == NULL)
                {
                  perror ("Allocating block memory");
                  exit (-1);
                }


              //  Start all "num_threads" threads to compute the mask.

              uint8_t *complete = (uint8_t *) calloc (num_threads, sizeof (uint8_t));
              if (complete == NULL)
                {
                  perror ("Allocating complete memory");
                  exit (-1);
                }

              double slat = (double) lat;
              double slon = (double) lon;

              for (int32_t i = 0 ; i < num_threads ; i++)
                {
                  mask_thread[i].mask (block, resolution, num_poly, poly_count, poly_y, poly_x, slat, slon, complete, num_threads, i);
                }


              //  We can't move on until all of the threads are complete.

              for (int32_t i = 0 ; i < num_threads ; i++)
                {
                  mask_thread[i].wait ();
                }


              int32_t size = (point_count * point_count) / 8;
              if ((point_count * point_count) % 8) size++;


              bit_block = (uint8_t *) calloc (size, sizeof (uint8_t));

              if (bit_block == NULL)
                {
                  perror ("Allocating bit_block memory");
                  exit (-1);
                }


              //  Copy the uint8_t block to the bit_block.

              int32_t block_size = point_count * point_count;
              for (int32_t pos = 0 ; pos < block_size ; pos++)
                {
                  if (block[pos])
                    {
                      bit_pack (bit_block, pos, 1, 1);
                    }
                  else
                    {
                      bit_pack (bit_block, pos, 1, 0);
                    }
                }


              //  Compress using zlib.

              uLong in_size = size;
              uLongf out_size = in_size + in_size * 0.10 + 100;
              uint8_t *out_buf = (uint8_t *) calloc (out_size, sizeof (uint8_t));

              if (out_buf == NULL)
                {
                  perror ("Allocating out_buf");
                  exit (-1);
                }


              int32_t n = compress2 (out_buf, &out_size, bit_block, size, 9);
              if (n)
                {
                  fprintf (stderr, "Error %d compressing record\n", n);
                  exit (-1);
                }


              //  Write the block address and the block size to the map.

              fseek (ofp, 0, SEEK_END);
              int32_t block_address = ftell (ofp);
              bit_pack (mapbuf, 0, 32, block_address);
              bit_pack (mapbuf, 32, 24, out_size);
              fseek (ofp, map_address, SEEK_SET);
              fwrite (mapbuf, 7, 1, ofp);

              if (lat == 10 && lon == 124) fprintf (stderr,"%s %s %d %d %ld %ld\n",NVFFL,block_address,out_size,in_size);

              //  Now write the block.

              fseek (ofp, block_address, SEEK_SET);
              fwrite (out_buf, out_size, 1, ofp);


              total_blocks++;
              total_block_size += (double) out_size;
              double avg = total_block_size / (double) total_blocks;

              fprintf (stderr, "%d blocks, average block size = %.2f\n\n", total_blocks, avg);


              /*
              //  Uncompressing test.

              in_size = (point_count * point_count) / 8 + 2000;

              fseek (ofp, block_address, SEEK_SET);
              fread (out_buf, out_size, 1, ofp);

              uint8_t *bit_box = (uint8_t *) calloc (in_size, sizeof (uint8_t));
              if (bit_box == NULL)
                {
                  perror ("Allocating bit_box memory in read_swbd_mask");
                  exit (-1);
                }

              int32_t status = uncompress (bit_box, &in_size, out_buf, out_size);
              if (status)
                {
                  fprintf (stderr, "Error %d uncompressing record\n", status);
                  exit (-1);
                }


              //  Copy the bit_block to the uint8_t block.

              for (int32_t pos = 0 ; pos < block_size ; pos++)
                {
                  block[pos] = (uint8_t) bit_unpack (bit_box, pos, 1);
                }

              free (bit_box);
              */


              free (out_buf);
              free (bit_block);
              free (block);


              //  Free all of the polygon memory.

              for (int32_t i = 0 ; i < num_poly ; i++)
                {
                  if (poly_x[i] != NULL) free (poly_x[i]);
                  if (poly_y[i] != NULL) free (poly_y[i]);
                }


              if (poly_x != NULL) free (poly_x);
              if (poly_y != NULL) free (poly_y);
              if (poly_count != NULL) free (poly_count);
              poly_x = NULL;
              poly_y = NULL;
              poly_count = NULL;
            }
        }
      //}
    }


  fclose (ofp);


  fprintf (stderr, "100%% processed                         \n\n");
  fflush (stderr);

  return (0);
}
