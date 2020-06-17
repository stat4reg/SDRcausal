//     ------------------------------------------------------------------------
//
//     This file is part of SDRcausal.
//
//     SDRcausal is free software: you can redistribute it and/or modify it
//     under the terms of the GNU General Public License as published by the
//     Free Software Foundation, either version 3 of the License, or (at your
//     option) any later version.
//
//     SDRcausal is distributed in the hope that it will be useful, but WITHOUT
//     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//     for more details.
//
//     You should have received a copy of the GNU General Public License along
//     with SDRcausal.  If not, see <https://www.gnu.org/licenses/>.
//
//     ------------------------------------------------------------------------
//
#include "index_element.h"

/** @brief Function used by qsort. See qsort manual for more information. */
int
compare_IndexElement (const void *a,
                      const void *b)
{
  if (((IndexElement*)a)->value > ((IndexElement*)b)->value)
    return 1;
  else if (((IndexElement*)a)->value < ((IndexElement*)b)->value)
    return -1;
  else
    return 0;
}
