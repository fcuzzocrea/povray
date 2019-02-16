//******************************************************************************
///
/// @file core/math/polynomialsolver.cpp
///
/// Implementations related to solving polynomial equations.
///
/// @author Alexander Enzmann
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

// Unit header file must be the first file included within POV-Ray *.cpp files (pulls in config)
#include "core/math/polynomialsolver.h"

#include "core/support/statistics.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

/*****************************************************************************
* Local preprocessor defines
******************************************************************************/

/// @var kSolveQuarticV1_Factor2
/// @brief const DBL value defining how close quartic equation is to being a square
/// of a quadratic.
///
/// @note
///     The closer this can get to zero before roots disappear, the less the chance
///     you will get spurious roots.
///
/// @attention
///     Used only in the old unused first version of solve_quartic().
///
const DBL kSolveQuarticV1_Factor2 = -1.0e-5;

/// @var kSolveQuarticV1_Factor3
/// @brief const DBL value similar to @ref kSolveQuarticV1_Factor2 at a later
/// stage of the algebraic solver.
///
/// @ref kSolveQuarticV1_Factor2 and @ref kSolveQuarticV1_Factor3 have been
/// defined so that quartic equations will properly render on fpu/compiler
/// combinations that only have 64 bit IEEE precision. Do not arbitrarily change
/// any of these values.
///
/// If you have a machine with a proper fpu/compiler combo - like a Mac with a
/// 68881, then use the native floating format (96 bits) and tune the values for
/// a little less tolerance: something like: factor2 = -1.0e-7, factor3 =
/// 1.0e-10. Twenty five years later the reality is still double accuracy
/// due use of fastmath (not IEEE compliant) compiling, use of SSE Fused
/// Multiply Add instructions, etc.
///
/// @attention
///     Used only in the old unused first version of solve_quartic().
///
const DBL kSolveQuarticV1_Factor3 = 1.0e-7;

/// @var kSolveQuarticV2_Factor4
/// @brief const DBL value used in the active solve_quartic() function.
///
/// Roughly acts as @ref kSolveQuarticV1_Factor2 and @ref
/// kSolveQuarticV1_Factor3 did for the original solve_quartic() versions,
/// but for the current solve_quartic() code.
///
const DBL kSolveQuarticV2_Factor4 = 1.0e-5;

/// @var kSolveCubic_2MultPiDiv3
/// const DBL value used in solve_cubic() equal to 2.0 * pi / 3.0.
///
const DBL kSolveCubic_2MultPiDiv3  = 2.0943951023931954923084;

/// @var kSolveCubic_4MultPiDiv3
/// const DBL value used in solve_cubic() equal to 4.0 * pi / 3.0.
///
const DBL kSolveCubic_4MultPiDiv3 = 4.1887902047863909846168;

/// @var kMaxIterations
/// const int max number of polysolve sturm chain based bisections.
///
/// @note
///     regula_falsa() uses twice this value internally as it can be
///     quite slow to converge in the worst case.
///
const int kMaxIterations = 65;

/// @var kSbisectMultRootThreshold
/// const @ref PRECISE_FLOAT value below which multiple roots ignored in sturm
/// chained based bisection and a single root at the middle of the current
/// interval is returned.
///
/// @note
///     Rays near tangent to surface create extremely close roots and instability in
///     sturm chain sign change results from numchanges(). Threshold often tripped
///     where the roots collapse inward due external equation set up creating a
///     multiple root at one t. The lathe and sphere_sweep objects being common
///     "offenders." Essentially, the polynomial curve is just touching the value=0
///     in these cases and not crossing through so as to create root 'multiplicity'
///     no root isolation method handles perfectly and with which solvers struggle.
///     The CEY early exits added to avoid a 'bug' are really this numerical issue
///     and yes, the exit provide by this value avoids a class of those early
///     returns with no roots due the CEY 'patch.' Might eventually consider some
///     other action, but that patch left as is for now.
///
const PRECISE_FLOAT kSbisectMultRootThreshold = (PRECISE_FLOAT)1e-10;

/// @var kRegulaFalsaThreshold
/// const @ref PRECISE_FLOAT threshold below which regula_falsa function is tried
/// when there is a single root.
///
/// @note
///     Ray interval max_value - min_value threshold below which regula_falsa
///     function is tried when there is a single root. Single roots happen often.
///     Rays continued by transparency or internal reflection for example will have
///     just one root. Idea is to delay use of regula_falsa method until the ray
///     domain is small given regula-falsi method can converge very, very slowly
///     with common enough ray-surface equations.
///
/// @todo
///     Initial setting 1.0 can likely be tuned for better performance trading off
///     bisection against regula_falsa.
///
const PRECISE_FLOAT kRegulaFalsaThreshold = (PRECISE_FLOAT)1.0;

/// @var kRelativeError
/// const @ref PRECISE_FLOAT smallest relative error along the ray when using
/// the polysolve(), sturm chain bisection / regula-falsi method.
///
const PRECISE_FLOAT kRelativeError = (PRECISE_FLOAT)1.0e-12;

/// @var kSolveQuadratic_SmallEnough
/// const @ref DBL value used to filter determinant value in older
/// solve_quadratic() in an unusual way causing artifacts.
///
const DBL kSolveQuadratic_SmallEnough = 1.0e-10;

///@note
///    The Newton-Raphson root polishing code herein is set up to loop to better solution.
///    However, the test criterea almost always does just one polishing pass where @ref DBL
///    set to double. Perhaps later refine the code.

/*****************************************************************************
* Local typedefs
******************************************************************************/

struct polynomial
{
    int ord;
    PRECISE_FLOAT coef[MAX_ORDER+1];
};


/*****************************************************************************
* Static functions
******************************************************************************/

static int modp (const polynomial *u, const polynomial *v, polynomial *r,
                 PRECISE_FLOAT zeroThreshold);
static bool regula_falsa (const int order, const PRECISE_FLOAT *coef, PRECISE_FLOAT a,
                          PRECISE_FLOAT b, DBL *val);
static int sbisect (int np, const polynomial *sseq, PRECISE_FLOAT min, PRECISE_FLOAT max,
                    int atmin, int atmax, DBL *roots);
static int numchanges (int np, const polynomial *sseq, PRECISE_FLOAT a);
static PRECISE_FLOAT polyeval (PRECISE_FLOAT x, int n, const PRECISE_FLOAT *Coeffs);
static int buildsturm (int ord, polynomial *sseq, PRECISE_FLOAT zeroThreshold);
static int visible_roots (int np, const polynomial *sseq);


/*****************************************************************************
*
* FUNCTION
*
*   modp
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Calculates the modulus of u(x) / v(x) leaving it in r.
*   It returns 0 if r(x) is a constant.
*
*   NOTE: This function assumes the leading coefficient of v is 1 or -1.
*         In more common terms that v is a monic polynomial.
*
* CHANGES
*
*   -
*
******************************************************************************/

static int modp(const polynomial *u, const polynomial *v, polynomial *r,
                PRECISE_FLOAT zeroThreshold)
{
    int k, j;

    *r = *u;

    if (v->coef[v->ord] < (PRECISE_FLOAT)0.0)
    {
        for (k = u->ord - v->ord - 1; k >= 0; k -= 2)
        {
            r->coef[k] = -r->coef[k];
        }

        for (k = u->ord - v->ord; k >= 0; k--)
        {
            for (j = v->ord + k - 1; j >= k; j--)
            {
                r->coef[j] = -r->coef[j] - r->coef[v->ord + k] * v->coef[j - k];
            }
        }
    }
    else
    {
        for (k = u->ord - v->ord; k >= 0; k--)
        {
            for (j = v->ord + k - 1; j >= k; j--)
            {
                r->coef[j] -= r->coef[v->ord + k] * v->coef[j - k];
            }
        }
    }

    k = v->ord - 1;

    // NOTE: Want a default value close to minimum 'effective' zeros as calculated
    // internally. Incoming 'effective' zeros should already have been set to
    // exactly 0.0. If an incoming polynomial is well behaved - no near zero values
    // during polynomial division - much smaller values can be used, but knowing
    // this a priori impossible. Too small a value can lead to false roots.
    //
    while (k >= 0 &&
           PRECISE_FABS(r->coef[k]) < zeroThreshold)
    {
        r->coef[k] = (PRECISE_FLOAT)0.0;

        k--;
    }

    r->ord = (k < 0) ? 0 : k;

    return(r->ord);
}



/*****************************************************************************
*
* FUNCTION
*
*   buildsturm
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Build the sturmian sequence for a polynomial.
*
* CHANGES
*
*   -
*
******************************************************************************/

static int buildsturm(int ord, polynomial *sseq, PRECISE_FLOAT zeroThreshold)
{
    int i;
    PRECISE_FLOAT f;
    PRECISE_FLOAT *fp, *fc;
    polynomial *sp;

    sseq[0].ord = ord;
    sseq[1].ord = ord - 1;

    /* calculate the derivative and normalize the leading coefficient. */

    f = PRECISE_FABS(sseq[0].coef[ord] * ord);

    fp = sseq[1].coef;
    fc = sseq[0].coef + 1;

    for (i = 1; i <= ord; i++)
    {
        *fp++ = *fc++ * i / f;
    }

    /* construct the rest of the Sturm sequence */

    for (sp = sseq + 2; modp(sp - 2, sp - 1, sp, zeroThreshold); sp++)
    {
        /* reverse the sign and normalize */

        f = -PRECISE_FABS(sp->coef[sp->ord]);

        for (fp = &sp->coef[sp->ord]; fp >= sp->coef; fp--)
        {
            *fp /= f;
        }
    }

    /* reverse the sign */

    sp->coef[0] = -sp->coef[0];

    return(sp - sseq);
}



/*****************************************************************************
*
* FUNCTION
*
*   visible_roots
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Find out how many visible intersections there are.
*
* CHANGES
*
*   -
*
******************************************************************************/

static int visible_roots(int np, const polynomial *sseq)
{
    int atposinf, atzero;
    const polynomial *s;
    PRECISE_FLOAT f, lf;

    atposinf = atzero = 0;

    /* changes at positve infinity */

    lf = sseq[0].coef[sseq[0].ord];

    for (s = sseq + 1; s <= sseq + np; s++)
    {
        f = s->coef[s->ord];

        if (lf == (PRECISE_FLOAT)0.0 || lf * f < (PRECISE_FLOAT)0.0)
        {
            atposinf++;
        }

        lf = f;
    }

    /* Changes at zero */

    lf = sseq[0].coef[0];

    for (s = sseq+1; s <= sseq + np; s++)
    {
        f = s->coef[0];

        if (lf == (PRECISE_FLOAT)0.0 || lf * f < (PRECISE_FLOAT)0.0)
        {
            atzero++;
        }

        lf = f;
    }

    return(atzero - atposinf);
}



/*****************************************************************************
*
* FUNCTION
*
*   numchanges
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Return the number of sign changes in the Sturm sequence in
*   sseq at the value a.
*
* CHANGES
*
*   -
*
******************************************************************************/

static int numchanges(int np, const polynomial *sseq, PRECISE_FLOAT a)
{
    int changes = 0;
    PRECISE_FLOAT f, lf, lfNotZero;
    const polynomial *s;

    lf = polyeval(a, sseq[0].ord, sseq[0].coef);
    lfNotZero = 0.0;

    for (s = sseq + 1; s <= sseq + np; s++)
    {
        f = polyeval(a, s->ord, s->coef);

        // Original code run 25 years counted multiple zeros as a sign change.
        // Almost never happens it seems, but zeros should just be skipped and
        // sign change check made to last non zero value seen.
        if (lf != (PRECISE_FLOAT)0.0)
        {
            if (lf * f < (PRECISE_FLOAT)0.0)
            {
                changes++;
            }
            lfNotZero = lf;
        }
        else
        {
            if (lfNotZero * f < (PRECISE_FLOAT)0.0)
            {
                changes++;
            }
        }
        lf = f;
    }

    return(changes);
}



/*****************************************************************************
*
* FUNCTION
*
*   sbisect
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Uses a bisection based on the sturm sequence for the polynomial
*   described in sseq to isolate intervals in which roots occur,
*   the roots are returned in the roots array in order of magnitude.
*
*   NOTE: This routine has one severe bug: When the interval containing the
*         root [min, max] has a root at one of its endpoints, as well as one
*         within the interval, the root at the endpoint will be returned
*         rather than the one inside. October, 2018. Saw no indication of
*         exactly this due sbisect though the regula_falsa code prior to recent
*         changes often returned false roots at ends over real roots if
*         evaluated values very small. The sturm chain method itself sees
*         roots on intervals of (minv, maxv] to older literature and testing.
*         This open low side has the potential to cause issues not in today's
*         default 0 to upper bound starting case, but should we support a
*         lower starting bound other than zero(1). In addition low side easy
*         'fixes' with reasonable performance would themselves be prone to
*         creating false roots due being value based.
*
*  (1) - Such capability would turn our Sturm chain based solver into a
*  general polynomial solver vs. being a specialized ray tracing one. Further,
*  min bound, not zero, methods are available which might lead to better
*  performance if we could pass a lower bound value >0.
*
* CHANGES
*
*   -
*
******************************************************************************/

static int sbisect(int np, const polynomial *sseq, PRECISE_FLOAT min_value, PRECISE_FLOAT max_value, int atmin, int atmax, DBL *roots)
{
    PRECISE_FLOAT mid;
    int n1, n2, its, atmid;

    if (((atmin - atmax) == 1) && ((max_value - min_value) < kRegulaFalsaThreshold))
    {
        /* first try using regula-falsa to find the root.  */

        if (regula_falsa(sseq->ord, sseq->coef, min_value, max_value, roots))
        {
            return(1);
        }
        else
        {
            /* That failed, so now find it by bisection */

            for (its = 0; its < kMaxIterations; its++)
            {
                mid = (min_value + max_value) / (PRECISE_FLOAT)2.0;

                atmid = numchanges(np, sseq, mid);

                /* The follow only happens if there is a bug.  And
                   unfortunately, there is. CEY 04/97
                 */
                if ((atmid<atmax) || (atmid>atmin))
                {
                    return(0);
                }

                if (PRECISE_FABS(mid) > kRelativeError)
                {
                    if (PRECISE_FABS((max_value - min_value) / mid) < kRelativeError)
                    {
                        roots[0] = (DBL)mid;

                        return(1);
                    }
                }
                else
                {
                    if (PRECISE_FABS(max_value - min_value) < kRelativeError)
                    {
                        roots[0] = (DBL)mid;

                        return(1);
                    }
                }

                if ((atmin - atmid) == 0)
                {
                    min_value = mid;
                }
                else
                {
                    max_value = mid;
                }
            }

            /* Bisection took too long - just return what we got */

            roots[0] = (DBL)mid;

            return(1);
        }
    }

    /* There is more than one root in the interval.
       Bisect to find new intervals. */

    for (its = 0; its < kMaxIterations; its++)
    {
        mid = (min_value + max_value) / (PRECISE_FLOAT)2.0;

        atmid = numchanges(np, sseq, mid);

        /* The follow only happens if there is a bug.  And
           unfortunately, there is. CEY 04/97
         */
        if ((atmid<atmax) || (atmid>atmin))
        {
            return(0);
        }

        if (PRECISE_FABS(mid) > kRelativeError)
        {
            if (PRECISE_FABS((max_value - min_value) / mid) < kRelativeError)
            {
                roots[0] = (DBL)mid;

                return(1);
            }
        }
        else
        {
            if (PRECISE_FABS(max_value - min_value) < kRelativeError)
            {
                roots[0] = (DBL)mid;

                return(1);
            }
        }

        n1 = atmin - atmid;
        n2 = atmid - atmax;

        if ((n1 != 0) && (n2 != 0))
        {
            if ((max_value - min_value) < kSbisectMultRootThreshold)
            {
                roots[0] = (DBL)mid;
                return(1);
            }

            n1 = sbisect(np, sseq, min_value, mid, atmin, atmid, roots);
            n2 = sbisect(np, sseq, mid, max_value, atmid, atmax, &roots[n1]);

            return(n1 + n2);
        }
        else
        {
            if ((n1 == 1) && (n2 == 0) && ((mid - min_value) < kRegulaFalsaThreshold))
            {
                n1 = sbisect(np, sseq, min_value, mid, atmin, atmid, roots);

                return(n1);
            }
            if ((n1 == 0) && (n2 == 1) && ((max_value - mid) < kRegulaFalsaThreshold))
            {
                n2 = sbisect(np, sseq, mid, max_value, atmid, atmax, roots);

                return(n2);
            }
        }

        if (n1 == 0)
        {
            min_value = mid;
        }
        else
        {
            max_value = mid;
        }
    }

    /* Took too long to bisect.  Just return what we got. */

    roots[0] = (DBL)mid;

    return(1);
}


/*****************************************************************************
*
* FUNCTION
*
*   polyeval
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Evaluate the value of a polynomial at the given value x.
*
*   The coefficients are stored in c in the following order:
*
*     c[0] + c[1] * x + c[2] * x ^ 2 + c[3] * x ^ 3 + ...
*
* CHANGES
*
*   -
*
******************************************************************************/

static PRECISE_FLOAT polyeval(PRECISE_FLOAT x, int n, const PRECISE_FLOAT *Coeffs)
{
    int i;
    PRECISE_FLOAT val;

    val = Coeffs[n];

    for (i = n-1; i >= 0; i--)
    {
        val = val * x + Coeffs[i];
    }

    return(val);
}



/*****************************************************************************
*
* FUNCTION
*
*   regular_falsa
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Close in on a root by using regula-falsa.
*
* CHANGES
*
*   -
*
******************************************************************************/

static bool regula_falsa(const int order, const PRECISE_FLOAT *coef, PRECISE_FLOAT a, PRECISE_FLOAT b, DBL *val)
{
    bool found=false;
    int its;
    PRECISE_FLOAT fa, fb, x, fx, lfx, mid;

    fa = polyeval(a, order, coef);
    fb = polyeval(b, order, coef);

    if (fa * fb > (PRECISE_FLOAT)0.0)
    {
        return(false);
    }

    lfx = fa;

    mid = (a + b) / (PRECISE_FLOAT)2.0;

    // NOTE: 2x kMaxIterations multiplier over sturm bisection requirement found to
    // happen in practice. Necessary if you want regula_falsa to find the root
    // where the algorithm converges very slowly.

    for (its = 0; its < (kMaxIterations * 2); its++)
    {
        x = (fb * a - fa * b) / (fb - fa);

        fx = polyeval(x, order, coef);

        if (PRECISE_FABS(mid) > kRelativeError)
        {
            if (PRECISE_FABS((b - a) / mid) < kRelativeError)
            {
                *val = (DBL)mid;

                found = true;

                break;
            }
        }
        else
        {
            if (PRECISE_FABS(b - a) < kRelativeError)
            {
                *val = (DBL)mid;

                found = true;

                break;
            }
            else if (PRECISE_FABS(fx) < gkPrecise_epsilon)
            {
                *val = (DBL)x;

                found = true;

                break;
            }
        }

        if (fa < (PRECISE_FLOAT)0.0)
        {
            if (fx < (PRECISE_FLOAT)0.0)
            {
                a = x;

                mid = (a + b) / (PRECISE_FLOAT)2.0;

                fa = fx;

                if ((lfx * fx) > (PRECISE_FLOAT)0.0)
                {
                    fb /= (PRECISE_FLOAT)2.0;
                }
            }
            else
            {
                b = x;

                mid = (a + b) / (PRECISE_FLOAT)2.0;

                fb = fx;

                if ((lfx * fx) > (PRECISE_FLOAT)0.0)
                {
                    fa /= (PRECISE_FLOAT)2.0;
                }
            }
        }
        else
        {
            if (fx < (PRECISE_FLOAT)0.0)
            {
                b = x;

                mid = (a + b) / (PRECISE_FLOAT)2.0;

                fb = fx;

                if ((lfx * fx) > (PRECISE_FLOAT)0.0)
                {
                    fa /= (PRECISE_FLOAT)2.0;
                }
            }
            else
            {
                a = x;

                mid = (a + b) / (PRECISE_FLOAT)2.0;

                fa = fx;

                if ((lfx * fx) > (PRECISE_FLOAT)0.0)
                {
                    fb /= (PRECISE_FLOAT)2.0;
                }
            }
        }

        lfx = fx;
    }

    return(found);
}



/*****************************************************************************
*
* FUNCTION
*
*   solve_quadratic
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Solve the quadratic equation:
*
*     x[0] * x^2 + x[1] * x + x[2] = 0.
*
*   The value returned by this function is the number of real roots.
*   The roots themselves are returned in y[0], y[1].
*
* CHANGES
*
*   -
*
******************************************************************************/

int solve_quadratic(const DBL *x, DBL *y)
{
#if (1)

    // NOTE: Numerically more accurate quadratic solver in a form from CERN
    // December 2016 presenation which is computationally more effecient
    // for most compilers - though gains differ substantially.
    //
    // See: "Portable SIMD and the VecCore Library." CERN. December 2016
    //
    // Numerical issues see:
    //     Numerical Recipes in C. 2nd Ed. or newer.
    //     https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
    //

    PRECISE_FLOAT a = x[0];
    PRECISE_FLOAT b = x[1];
    PRECISE_FLOAT c = x[2];
    PRECISE_FLOAT tmp;

    if (PRECISE_FABS(a) < gkDBL_epsilon)
    {
        if (PRECISE_FABS(b) < gkDBL_epsilon)
        {
            return(0);
        }

        y[0] = (DBL)(-c / b);

        return(1);
    }

    PRECISE_FLOAT a_inv = (PRECISE_FLOAT)1.0 / a;
    PRECISE_FLOAT delta = b * b - (PRECISE_FLOAT)4.0 * a * c;
    PRECISE_FLOAT s     = (b >= 0) ? (PRECISE_FLOAT)1.0 : (PRECISE_FLOAT)-1.0;

    int roots = delta > (PRECISE_FLOAT)0.0 ? 2 : delta < (PRECISE_FLOAT)0.0 ? 0 : 1;

    switch (roots) {
    case 2:
      tmp  = (PRECISE_FLOAT)-0.5 * (b + s * PRECISE_SQRT(delta));
      y[1] = (DBL)(c / tmp);
      y[0] = (DBL)(tmp * a_inv);
      return(roots);

    case 0:
      return(roots);

    case 1:
      y[0] = (DBL)((PRECISE_FLOAT)-0.5 * b * a_inv);
      return(roots);

    default:
      return(0);
    }

#else

    // NOTE: POV-Ray versions of the quadratic solver have been atypical since
    // the original 'traditional' implementation. The version below only works
    // reliably for roots where a and b values limited to >= about 1e-10.

    DBL d, t, a, b, c;

    a = x[0];
    b = -x[1];
    c = x[2];

    if (fabs(a) < kSolveQuadratic_SmallEnough)
    {
        if (fabs(b) < kSolveQuadratic_SmallEnough)
        {
            return(0);
        }

        y[0] = c / b;

        return(1);
    }

    // Normalize the coefficients. Added for v3.7.
    b /= a;
    c /= a;
    a  = 1.0;

    d = b * b - 4.0 * a * c;

    // Treat values of d around 0 as 0. Effectively returns one root
    // where supposing infinite numerical accuracy there are two.
    // The cost is bad roots - fuzz / inaccuracy in some results.

    if ((d > -kSolveQuadratic_SmallEnough) && (d < kSolveQuadratic_SmallEnough))
    {
        y[0] = 0.5 * b / a;

        return(1);
    }
    else
    {
        if (d < 0.0)
        {
            return(0);
        }
    }

    d = sqrt(d);

    t = 2.0 * a;

    y[0] = (b + d) / t;
    y[1] = (b - d) / t;

    return(2);
#endif
}



/*****************************************************************************
*
* FUNCTION
*
*   solve_cubic
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Solve the cubic equation:
*
*     x[0] * x^3 + x[1] * x^2 + x[2] * x + x[3] = 0.
*
*   The result of this function is an integer that tells how many real
*   roots exist.  Determination of how many are distinct is up to the
*   process that calls this routine.  The roots that exist are stored
*   in (y[0], y[1], y[2]).
*
*   NOTE: This function relies very heavily on trigonometric functions and
*         the square root function.  If an alternative solution is found
*         that does not rely on transcendentals this code will be replaced.
*
* CHANGES
*
*   -
*
******************************************************************************/

int solve_cubic(const DBL *x, DBL *y)
{
    DBL Q, R, Q3, R2, sQ, d, an, theta;
    DBL A1sqr, a0, a1, a2, a3;
    int rootCount;

    // @todo
    //     See GNU GSL library for a version which handles both root multiplicities
    //     and better real roots near the x axis. Mostly concerned because this
    //     solver used as part of solve_quartic where such 'clarity' of result
    //     may matter to that solver's result. What we have here - aside from accuracy
    //     problems requiring the root polishing - seems to work OK for ray tracing.
    //     It might be solve_quartic should use a gsl like solve_cubic implementation.
    //

    a0 = x[0];

    if (fabs(a0) < gkDBL_epsilon)
    {
        return(solve_quadratic(&x[1], y));
    }
    else
    {
        // Note. Code since v1.0 had conditional for a0 == 1 which would avoid the
        // division by a0 if possible. If in an environment where polynomials more
        // often in monic form this 'might' make sense. In POV-Ray where very often
        // not monic and today with SIMD hardware common, it does not.
        a1 = x[1] / a0;
        a2 = x[2] / a0;
        a3 = x[3] / a0;
    }

    A1sqr = a1 * a1;

    Q = (A1sqr - 3.0 * a2) / 9.0;

    // @todo
    //     Look to remove the DJGPP bug fix below. While on the surface it removes a couple
    //     muliplications it forfeits accuracy in the subtraction step within the parenthesis.

    /* Modified to save some multiplications and to avoid a floating point
       exception that occured with DJGPP and full optimization. [DB 8/94] */

    R = (a1 * (A1sqr - 4.5 * a2) + 13.5 * a3) / 27.0;

    Q3 = Q * Q * Q;

    R2 = R * R;

    d = Q3 - R2;

    an = a1 / 3.0;

    if (d >= 0.0)
    {
        /* Three real roots. */

        d = R / sqrt(Q3);

        theta = acos(d) / 3.0;

        sQ = -2.0 * sqrt(Q);

        y[0] = sQ * cos(theta) - an;
        y[1] = sQ * cos(theta + kSolveCubic_2MultPiDiv3) - an;
        y[2] = sQ * cos(theta + kSolveCubic_4MultPiDiv3) - an;

        rootCount = 3;
    }
    else
    {
        sQ = pow(sqrt(R2 - Q3) + fabs(R), 1.0 / 3.0);

        if (R < 0)
        {
            y[0] = (sQ + Q / sQ) - an;
        }
        else
        {
            y[0] = -(sQ + Q / sQ) - an;
        }

        rootCount = 1;
    }

    // Newton-Raphson root polishing.
    PRECISE_FLOAT pv, pvLast, dpv, dt, t;

    auto lambdaCubicPolyValAnDerVal = [&pv,&dpv](const PRECISE_FLOAT t, const DBL *x) -> void
    {
         pv  = (PRECISE_FLOAT)x[0] * t + (PRECISE_FLOAT)x[1];
         dpv = (PRECISE_FLOAT)x[0];
         for (int k=2; k<=3; k++)
         {
             dpv = dpv * t + pv;
             pv  = pv * t + (PRECISE_FLOAT)x[k];
         }
    };

    for (int c = 0; c < rootCount; c++)
    {
        if (y[c] <= 0.0)
        {
            continue;
        }
        t = (PRECISE_FLOAT)y[c];
        lambdaCubicPolyValAnDerVal(t,&x[0]);
        pvLast = pv;

        for (int j = 0; j < 7; j++)
        {
            if (dpv == 0.0)
            {
               break;
            }
            else if (pv == 0.0)
            {
                y[c] = (DBL)t;
                break;
            }
            dt = pv / dpv;
            t -= dt;

            lambdaCubicPolyValAnDerVal(t,&x[0]);

            if (PRECISE_FABS(dt) < gkMinIsectDepthReturned)
            {
                if (PRECISE_FABS(pv) < PRECISE_FABS(pvLast))
                {
                    y[c] = (DBL)t;
                }
                break;
            }
        }
    }
    return(rootCount);
}


#if (0)
/*****************************************************************************
*
* FUNCTION
*
*   solve_quartic
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   The old way of solving quartics algebraically.
*   This is an adaptation of the method of Lodovico Ferrari (Circa 1731).
*
* CHANGES
*
*   -
*
******************************************************************************/

int solve_quartic(const DBL *x, DBL *results)
{
    DBL cubic[4], roots[3];
    DBL a0, a1, y, d1, x1, t1, t2;
    DBL c0, c1, c2, c3, c4, d2, q1, q2;
    int i;

    c0 = x[0];

    if (fabs(c0) < gkDBL_epsilon)
    {
        return(solve_cubic(&x[1], results));
    }
    else
    {
        if (c0 != 1.0)
        {
            c1 = x[1] / c0;
            c2 = x[2] / c0;
            c3 = x[3] / c0;
            c4 = x[4] / c0;
        }
        else
        {
            c1 = x[1];
            c2 = x[2];
            c3 = x[3];
            c4 = x[4];
        }
    }

    /* The first step is to take the original equation:

         x^4 + b*x^3 + c*x^2 + d*x + e = 0

       and rewrite it as:

         x^4 + b*x^3 = -c*x^2 - d*x - e,

       adding (b*x/2)^2 + (x^2 + b*x/2)y + y^2/4 to each side gives a
       perfect square on the lhs:

         (x^2 + b*x/2 + y/2)^2 = (b^2/4 - c + y)x^2 + (b*y/2 - d)x + y^2/4 - e

       By choosing the appropriate value for y, the rhs can be made a perfect
       square also.  This value is found when the rhs is treated as a quadratic
       in x with the discriminant equal to 0.  This will be true when:

         (b*y/2 - d)^2 - 4.0 * (b^2/4 - c*y)*(y^2/4 - e) = 0, or

         y^3 - c*y^2 + (b*d - 4*e)*y - b^2*e + 4*c*e - d^2 = 0.

       This is called the resolvent of the quartic equation.  */

    a0 = 4.0 * c4;

    cubic[0] = 1.0;
    cubic[1] = -1.0 * c2;
    cubic[2] = c1 * c3 - a0;
    cubic[3] = a0 * c2 - c1 * c1 * c4 - c3 * c3;

    i = solve_cubic(&cubic[0], &roots[0]);

    if (i > 0)
    {
        y = roots[0];
    }
    else
    {
        return(0);
    }

    /* What we are left with is a quadratic squared on the lhs and a
       linear term on the right.  The linear term has one of two signs,
       take each and add it to the lhs.  The form of the quartic is now:

         a' = b^2/4 - c + y,    b' = b*y/2 - d, (from rhs quadritic above)

         (x^2 + b*x/2 + y/2) = +sqrt(a'*(x + 1/2 * b'/a')^2), and
         (x^2 + b*x/2 + y/2) = -sqrt(a'*(x + 1/2 * b'/a')^2).

       By taking the linear term from each of the right hand sides and
       adding to the appropriate part of the left hand side, two quadratic
       formulas are created.  By solving each of these the four roots of
       the quartic are determined.
    */

    i = 0;

    a0 = c1 / 2.0;
    a1 = y / 2.0;

    t1 = a0 * a0 - c2 + y;

    if (t1 < 0.0)
    {
        if (t1 > kSolveQuarticV1_Factor2)
        {
            t1 = 0.0;
        }
        else
        {
            /* First Special case, a' < 0 means all roots are complex. */

            return(0);
         }
     }

    if (t1 < kSolveQuarticV1_Factor3)
    {
        /* Second special case, the "x" term on the right hand side above
           has vanished.  In this case:

             (x^2 + b*x/2 + y/2) = +sqrt(y^2/4 - e), and
             (x^2 + b*x/2 + y/2) = -sqrt(y^2/4 - e).  */

        t2 = a1 * a1 - c4;

        if (t2 < 0.0)
        {
            return(0);
        }

        x1 = 0.0;
        d1 = sqrt(t2);
    }
    else
    {
        x1 = sqrt(t1);
        d1 = 0.5 * (a0 * y - c3) / x1;
    }

    /* Solve the first quadratic */

    q1 = -a0 - x1;
    q2 = a1 + d1;
    d2 = q1 * q1 - 4.0 * q2;

    if (d2 >= 0.0)
    {
        d2 = sqrt(d2);

        results[0] = 0.5 * (q1 + d2);
        results[1] = 0.5 * (q1 - d2);

        i = 2;
    }

    /* Solve the second quadratic */

    q1 = q1 + x1 + x1;
    q2 = a1 - d1;
    d2 = q1 * q1 - 4.0 * q2;

    if (d2 >= 0.0)
    {
        d2 = sqrt(d2);

        results[i++] = 0.5 * (q1 + d2);
        results[i++] = 0.5 * (q1 - d2);
    }

    return(i);
}
#else
/*****************************************************************************
*
* FUNCTION
*
*   solve_quartic
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Solve a quartic using the method of Francois Vieta (Circa 1735).
*
* CHANGES
*
*   -
*
******************************************************************************/

int solve_quartic(const DBL *x, DBL *results)
{
    DBL cubic[4];
    DBL roots[3];
    DBL c12, z, p, q, q1, q2, r, d1, d2;
    DBL c0, c1, c2, c3, c4;
    int i;

    /* Make sure the quartic has a leading coefficient of 1.0 */

    c0 = x[0];

    if (fabs(c0) < gkDBL_epsilon)
    {
        return(solve_cubic(&x[1], results));
    }
    else
    {
        // Note. Code since v1.0 had conditional for c0 == 1 which would avoid the
        // division by c0 if possible. If in an environment where polynomials more
        // often in monic form this 'might' make sense. In POV-Ray where very often
        // not monic and today with SIMD hardware common, it does not.
        c1 = x[1] / c0;
        c2 = x[2] / c0;
        c3 = x[3] / c0;
        c4 = x[4] / c0;
    }

    /* Compute the cubic resolvant */

    c12 = c1 * c1;
    p = -0.375 * c12 + c2;
    q = 0.125 * c12 * c1 - 0.5 * c1 * c2 + c3;
    r = -0.01171875 * c12 * c12 + 0.0625 * c12 * c2 - 0.25 * c1 * c3 + c4;

    cubic[0] = 1.0;
    cubic[1] = -0.5 * p;
    cubic[2] = -r;
    cubic[3] = 0.5 * r * p - 0.125 * q * q;

    i = solve_cubic(cubic, roots);

    if (i > 0)
    {
        z = roots[0];
    }
    else
    {
        return(0);
    }

    d1 = 2.0 * z - p;

    if (d1 < 0.0)
    {
        if (d1 > -kSolveQuarticV2_Factor4)
        {
            d1 = 0.0;
        }
        else
        {
            return(0);
        }
    }

    if (d1 < kSolveQuarticV2_Factor4)
    {
        d2 = z * z - r;

        if (d2 < 0.0)
        {
            return(0);
        }

        d2 = sqrt(d2);
    }
    else
    {
        d1 = sqrt(d1);
        d2 = 0.5 * q / d1;
    }

    /* Set up useful values for the quadratic factors */

    q1 = d1 * d1;
    q2 = -0.25 * c1;

    i = 0;

    /* Solve the first quadratic */

    p = q1 - 4.0 * (z - d2);

    if (p == 0)
    {
        results[i++] = -0.5 * d1 - q2;
    }
    else
    {
        if (p > 0)
        {
            p = sqrt(p);
            results[i++] = -0.5 * (d1 + p) + q2;
            results[i++] = -0.5 * (d1 - p) + q2;
        }
    }

    /* Solve the second quadratic */

    p = q1 - 4.0 * (z + d2);

    if (p == 0)
    {
        results[i++] = 0.5 * d1 - q2;
    }
    else
    {
        if (p > 0)
        {
            p = sqrt(p);
            results[i++] = 0.5 * (d1 + p) + q2;
            results[i++] = 0.5 * (d1 - p) + q2;
        }
    }

    // The quartic root finding method above somewhat inaccurate and amplified
    // at large scale up. The follow code polishes the found roots so, for example,
    // <=0.0 roots which otherwise drift positive and cause artifacts are corrected.
    // This the reason for the historically too large 1e-2 minimum intersection
    // depth used with blobs that itself caused artifacts of other kinds. Further,
    // that blob 1e-2 value was not large enough to filter all drift to 0+ roots in
    // any case without the polishing below. Newton-Raphson method.
    //
    PRECISE_FLOAT pv, pvLast, dpv, dt, t;

    auto lambdaQuarticPolyValAnDerVal = [&pv,&dpv](const PRECISE_FLOAT t, const DBL *x) -> void
    {
         pv  = (PRECISE_FLOAT)x[0] * t + (PRECISE_FLOAT)x[1];
         dpv = (PRECISE_FLOAT)x[0];
         for (int k=2; k<=4; k++)
         {
             dpv = dpv * t + pv;
             pv  = pv * t + (PRECISE_FLOAT)x[k];
         }
    };

    for (int c = 0; c < i; c++)
    {
        if (results[c] <= 0.0)
        {
            continue;
        }
        t = (PRECISE_FLOAT)results[c];
        lambdaQuarticPolyValAnDerVal(t,&x[0]);
        pvLast = pv;

        for (int j = 0; j < 7; j++)
        {
            if (dpv == 0.0)
            {
               break;
            }
            else if (pv == 0.0)
            {
                results[c] = (DBL)t;
                break;
            }
            dt = pv / dpv;
            t -= dt;

            lambdaQuarticPolyValAnDerVal(t,&x[0]);

            if (PRECISE_FABS(dt) < gkMinIsectDepthReturned)
            {
                if (PRECISE_FABS(pv) < PRECISE_FABS(pvLast))
                {
                    results[c] = (DBL)t;
                }
                break;
            }
        }
    }

    return(i);
}
#endif



/*****************************************************************************
*
* FUNCTION
*
*   polysolve
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Root solver based on the Sturm sequences for a polynomial.
*
* CHANGES
*
*   -
*
******************************************************************************/

int polysolve (int order, const DBL *Coeffs, DBL *roots, DBL HandleCollapsedRootsValue, DBL MaxBound)
{
    polynomial sseq[MAX_ORDER+1];
    DBL min_value, max_value, max_value2, Abs_Coeff_n;
    int i, nroots, np, atmin, atmax;
    bool potentialLeadingZero = true;
    PRECISE_FLOAT PreciseFltValue;

    //---
    // Reverse coefficients into order used herein.
    // Drop any leading original order coefficients meaningfully 0.0.
    // (Commented below) Set any remaining coefficients meaningfully 0.0 to exactly 0.0.

    np = 0;
    for (i = 0; i <= order; i++)
    {
        if (potentialLeadingZero && fabs(Coeffs[i]) < gkDBL_epsilon)
        {
            np++;
        }
        else
        {
            potentialLeadingZero = false;
          //if (fabs(Coeffs[i]) < gkDBL_epsilon)
          //{
          //    sseq[0].coef[order-i] = (PRECISE_FLOAT)0.0;
          //}
          //else
          //{
                sseq[0].coef[order-i] = (PRECISE_FLOAT)Coeffs[i];
          //}
        }
    }
    order -= np;
    if (order == 0)
    {
        return(0);
    }
    else if (order == 1)
    {
        if (sseq[0].coef[1] != (PRECISE_FLOAT)0.0)
        {
            roots[0] = -sseq[0].coef[0] / sseq[0].coef[1];

            return(1);
        }
        else
        {
            return(0);
        }
    }

    //---
    // Build the Sturm sequence

    if (HandleCollapsedRootsValue > 0.0)
    {
        nroots = 0;

        PreciseFltValue = HandleCollapsedRootsValue;
        np = buildsturm(order, &sseq[0], PreciseFltValue);

        // At other than default as zero value to buildsturm above, visible_roots is sometimes
        // less accurate than numchanges. Using numchanges for nroots calc.
        atmin = numchanges(np, sseq, (PRECISE_FLOAT)0.0);
        if (MaxBound > 0.0)
        {
            max_value = max(1.0,MaxBound);
            atmax = numchanges(np, sseq, (PRECISE_FLOAT)max_value);
        }
        else
        {
            atmax = numchanges(np, sseq, (PRECISE_FLOAT)MAX_DISTANCE);
        }
        nroots = atmin - atmax;
    }
    else
    {
        // So long as incoming coefficients calculated at DBL, set the
        // efffective zero in following call to gkDBL_epsilon. If selected
        // shapes move to PRECISE_FLOAT, pass better HandleCollapsedRootsValue
        //
        np = buildsturm(order, &sseq[0], (PRECISE_FLOAT)gkDBL_epsilon);
        nroots = visible_roots(np, sseq);
    }

    // For total number of visible roots.
    // NOTE: Changed to <=0 test over ==0 due sphere_sweep b_spline
    // going negative when the modp leading coef filter set lower.
    // Similar change to the numchanges based test below.

    if (nroots <= 0)
    {
        return(0);
    }



    //---
    // Bracket the roots
    if (MaxBound > 0.0)
    {
        min_value = 0.0;
        max_value = max(1.0,MaxBound); // Prevent passed value being too small.

        atmin = numchanges(np, sseq, (PRECISE_FLOAT)min_value);
        atmax = numchanges(np, sseq, (PRECISE_FLOAT)max_value);
    }
    else
    {
        min_value = 0.0;

        // Use variation of Augustin-Louis Cauchy's method to determine an upper bound for max_value.

        // Tighter upper bound found at:
        //    https://en.wikipedia.org/wiki/Properties_of_polynomial_roots#Other_bounds
        // which took it from:
        //    Cohen, Alan M. (2009). "Bounds for the roots of polynomial equations".
        //    Mathematical Gazette. 93: 87-88.
        // NOTE: Had to use > 1.0 in max_value2 calculation in practice...

        Abs_Coeff_n = fabs(Coeffs[0]); // Leading zeros dropped above.
        max_value2  = 1.1 + fabs(Coeffs[1]/Abs_Coeff_n);
        max_value   = fabs(Coeffs[2]);
        for (i = 3; i <= order; i++)
        {
            max_value = max(fabs(Coeffs[i]),max_value);
        }
        max_value /= Abs_Coeff_n + EPSILON;
        max_value = min(max(max_value,max_value2),MAX_DISTANCE);

        // NOTE: Found in practice roots occasionally, slightly outside upper bound...
        // Perhaps related to how the sturm chain is pruned in modp(). Until sorted adding
        // the following sanity check which restores a MAX_DISTANCE upper bound where
        // root(s) exists above estimated upper bound.

        atmin = numchanges(np, sseq, (PRECISE_FLOAT)max_value);
        atmax = numchanges(np, sseq, (PRECISE_FLOAT)MAX_DISTANCE);
        if ((atmin - atmax) != 0)
        {
            max_value = MAX_DISTANCE;
        }
        else
        {
            atmax = atmin;
        }
        atmin = numchanges(np, sseq, (PRECISE_FLOAT)min_value);
    }

    nroots = atmin - atmax;

    if (nroots <= 0)
    {
        return(0);
    }

    // Perform the bisection.
    nroots = sbisect(np, sseq, (PRECISE_FLOAT)min_value, (PRECISE_FLOAT)max_value,
             atmin, atmax, roots);

    // Newton-Raphson root polishing step.
    PRECISE_FLOAT pv, pvLast, dpv, dt, t;

    auto lambdaPolysolvePolyValAnDerVal = [&pv,&dpv]
         (const int ord, const PRECISE_FLOAT t, const polynomial *sseq) -> void
    {
         pv  = (PRECISE_FLOAT)sseq[0].coef[ord] * t + (PRECISE_FLOAT)sseq[0].coef[ord-1];
         dpv = (PRECISE_FLOAT)sseq[0].coef[ord];
         for (int k=ord-2; k>=0; k--)
         {
             dpv = dpv * t + pv;
             pv  = pv * t + (PRECISE_FLOAT)sseq[0].coef[k];
         }
    };

    for (int c = 0; c < nroots; c++)
    {
        t = (PRECISE_FLOAT)roots[c];
        lambdaPolysolvePolyValAnDerVal(order,t,&sseq[0]);
        pvLast = pv;

        for (int j = 0; j < 7; j++)
        {
            if (dpv == 0.0)
            {
               break;
            }
            else if (pv == 0.0)
            {
                roots[c] = (DBL)t;
                break;
            }
            dt = pv / dpv;
            t -= dt;

            lambdaPolysolvePolyValAnDerVal(order,t,&sseq[0]);

            if (PRECISE_FABS(dt) < gkMinIsectDepthReturned)
            {
                if (PRECISE_FABS(pv) < PRECISE_FABS(pvLast))
                {
                    roots[c] = (DBL)t;
                }
                break;
            }
        }
    }

    return(nroots);
}

}
