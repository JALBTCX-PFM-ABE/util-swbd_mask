
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



#include "maskThread.hpp"

maskThread::maskThread (QObject *parent)
  : QThread(parent)
{
}



maskThread::~maskThread ()
{
}



void maskThread::mask (uint8_t *bl, int32_t r, int32_t np, int32_t *pc, double **py, double **px, double slt, double sln, uint8_t *c, int32_t nt,
                       int32_t p)
{
  QMutexLocker locker (&mutex);

  l_block = bl;
  l_resolution = r;
  l_num_poly = np;
  l_poly_count = pc;
  l_poly_y = py;
  l_poly_x = px;
  l_sw_lat = slt;
  l_sw_lon = sln;
  l_complete = c;
  l_num_threads = nt;
  l_pass = p;

  if (!isRunning ()) start ();
}



void maskThread::run ()
{
  mutex.lock ();

  uint8_t *block = l_block;
  int32_t resolution = l_resolution;
  int32_t num_poly = l_num_poly;
  int32_t *poly_count = l_poly_count;
  double **poly_y = l_poly_y;
  double **poly_x = l_poly_x;
  double sw_lat = l_sw_lat;
  double sw_lon = l_sw_lon;
  uint8_t *complete = l_complete;
  int32_t num_threads = l_num_threads;
  int32_t pass = l_pass;

  mutex.unlock ();


  //qDebug () << __LINE__ << pass;

  int32_t point_count = 3600 / resolution;
  int32_t percent = 0, old_percent = -1;
  int32_t block_count = NINT (sqrt ((double) num_threads));
  int32_t pass_point_count = point_count / block_count;

  int32_t start_x = (pass % block_count) * pass_point_count;
  int32_t start_y = (pass / block_count) * pass_point_count;

  /*
  switch (pass)
    {
    case 0:
      start_y = 0;
      start_x = 0;
      break;

    case 1:
      start_y = 0;
      start_x = pass_point_count;
      break;

    case 2:
      start_y = pass_point_count;
      start_x = 0;
      break;

    case 3:
      start_y = pass_point_count;
      start_x = pass_point_count;
      break;
    }
  */

  int32_t end_x = start_x + pass_point_count;
  int32_t end_y = start_y + pass_point_count;


  double pc_double = (double) point_count;
  double new_pc_double = (double) pass_point_count;


  //  Latitude loop.

  for (int32_t i = start_y ; i < end_y ; i++)
    {
      //  Compute the latitude of the center of the "spacing" sized bin (that's why we add 0.5).

      double slat = (double) sw_lat + (double) (i + 0.5) / pc_double;


      //  Longitude loop.

      for (int32_t j = start_x ; j < end_x ; j++)
        {
          //  Compute the longitude of the center of the "spacing" sized bin (that's why we add 0.5).

          double slon = (double) sw_lon + (double) (j + 0.5) / pc_double;


          int32_t inside_count = 0;


          //  Check against all polygons.

          for (int32_t k = 0 ; k < num_poly ; k++)
            {
              if (inside_polygon2 (poly_x[k], poly_y[k], poly_count[k], slon, slat)) inside_count++;
            }


          //  Set the flag (NVTrue for land, NVFalse for water).

          if (inside_count % 2)
            {
              block[i * point_count + j] = NVFalse;
            }
          else
            {
              block[i * point_count + j] = NVTrue;
            }
        }


      percent = (int32_t) (((double) (i - start_y) / new_pc_double) * 100.0);
      if (percent != old_percent)
        {
          fprintf (stderr, "Pass %d - %03d%% processed\r", pass, percent);
          fflush (stderr);
          old_percent = percent;
        }
    }


  //qDebug () << __LINE__ << pass;

  complete[pass] = NVTrue;
}
