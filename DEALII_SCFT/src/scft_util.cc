/*
 * scft_util.cc
 *
 *  Created on: Oct 8, 2018
 *      Author: chen
 */

#include "scft_urtil.h"
#include <stdio.h>
#include <stdlib.h>
#include<string>

void
read_yita_middle_1D (std::vector<double>& yita_middle, std::string filename,
		     int& N)
{

  FILE *file;
  file = fopen (filename.c_str (), "r");
  if (file == NULL)
    {
      fprintf (stderr, "Can't open input file!\n");
      exit (1);
    }
  char buff[255];
  int i;
  double value, x;
  if (fgets (buff, 255, file))
    {
      ;
    }
  sscanf (buff, "N= %d", &N);
  fgets (buff, 255, file);  // jump the line of mean_field_free_energy
  yita_middle.resize (N,0.);
  while (fgets (buff, 255, file))
    {
      sscanf (buff, "%d ,%lf, %lf", &i, &x, &value);
      yita_middle[i] = value;
    }
  fclose (file);
}

