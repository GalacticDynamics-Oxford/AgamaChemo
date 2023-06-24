/** \file    obs.h 
    \brief   Sky coordinates & functions relating to observations
    \author  Fred Thompson & James Binney
    \date    2021

The obs_ files provide code to facilitate interaction with observations.
   In obs_base polar sky coords (right acsension, declination) or Galactic (longitude, latitude)
   are defined alongside associated proper motions (with mul = dot\ell\cos(b), etc).
   A solarhifter moves data between Galactocentric and sky coordinates.
   The obs_dust files have code to describe dust distributions
   The obs_los files define lines of sight with functions to compute Galctocentric coords
   to/from distance down los. If a dustModel is given, xtinction along the los is computed 

*/

#pragma once
#include <cassert>
#include "coord.h"
#include "units.h"
#include "obs_base.h"
#include "obs_dust.h"
#include "obs_los.h"