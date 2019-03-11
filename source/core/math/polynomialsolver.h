//******************************************************************************
///
/// @file core/math/polynomialsolver.h
///
/// Declarations related to solving polynomial equations.
///
/// @copyright
/// @parblock
///
/// Persistence of Vision Ray Tracer ('POV-Ray') version 3.8.
/// Copyright 1991-2019 Persistence of Vision Raytracer Pty. Ltd.
///
/// POV-Ray is free software: you can redistribute it and/or modify
/// it under the terms of the GNU Affero General Public License as
/// published by the Free Software Foundation, either version 3 of the
/// License, or (at your option) any later version.
///
/// POV-Ray is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU Affero General Public License for more details.
///
/// You should have received a copy of the GNU Affero General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.
///
/// ----------------------------------------------------------------------------
///
/// POV-Ray is based on the popular DKB raytracer version 2.12.
/// DKBTrace was originally written by David K. Buck.
/// DKBTrace Ver 2.0-2.12 were written by David K. Buck & Aaron A. Collins.
///
/// @endparblock
///
//******************************************************************************

#ifndef POVRAY_CORE_POLYNOMIALSOLVER_H
#define POVRAY_CORE_POLYNOMIALSOLVER_H

// Module config header file must be the first file included within POV-Ray unit header files
#include "core/configcore.h"

// C++ variants of C standard header files
// C++ standard header files
//  (none at the moment)

// POV-Ray header files (base module)
// POV-Ray header files (core module)
//  (none at the moment)

namespace pov
{

//##############################################################################
///
/// @defgroup PovCoreMathPolynomialsolver Polynomial Solver
/// @ingroup PovCoreMath
///
/// @{

class RenderStatistics;

/*****************************************************************************
* Global preprocessor defines
******************************************************************************/

/// @def MAX_ORDER
/// Maximum supported polynomial order.
///
/// @todo
///     This currently carries a large, fixed, per polysolve() call memory allocation
///     on the stack. Size is on the order of (MAX_ORDER+1)*int + PRECISE_FLOAT *
///     (MAX_ORDER+1)^2 which impacts performance and performance stability
///     especially for threads > physical cores. Allocation based on current
///     equation order would be better.
///
#define MAX_ORDER 35

/// @todo
///     With the elimination of the Solve_Polynomial wrapper we are no longer
///     setting stats[Polynomials_Tested] and stats[Roots_Eliminated]. Note both
///     long inaccurate. The first never capturing anything close to the actual
///     polynomials solved for roots. Rather only calls to Solve_Polynomial. The
///     roots eliminated value was entirely inaccurate with respect to roots
///     eliminated and rather reflected one of a few ways coefficients were
///     being changed ahead of the actual solver work. Currently leaving the
///     stats infrastructure there with the thought both could be restored more
///     accurately, but truthfully doubt there would ever be any actual value to
///     users or developers.
///

/// @remark
///     The coefficient order for the polynomial when calling any function defined
///     herein is:
///
///         c[0] * x ^ n + c[1] * x ^ (n-1) + ... + c[n-1] * x + c[n] = 0
///

/*****************************************************************************
* Global functions
******************************************************************************/

int polysolve (int order, const DBL *Coeffs, DBL *roots,
               DBL HandleCollapsedRootsValue, DBL MaxBound);
int solve_quadratic (const DBL *x, DBL *y);
int solve_cubic (const DBL *x, DBL *y);
int solve_quartic (const DBL *x, DBL *y);

/// @}
///
//##############################################################################

}
// end of namespace pov

#endif // POVRAY_CORE_POLYNOMIALSOLVER_H
