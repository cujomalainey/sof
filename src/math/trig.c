// SPDX-License-Identifier: BSD-3-Clause
//
// Copyright(c) 2016 Intel Corporation. All rights reserved.
//
// Author: Seppo Ingalsuo <seppo.ingalsuo@linux.intel.com>
//         Liam Girdwood <liam.r.girdwood@linux.intel.com>
//         Keyon Jie <yang.jie@linux.intel.com>
//         Shriram Shastry <malladi.sastry@linux.intel.com>

#include <sof/audio/format.h>
#include <sof/math/trig.h>
#include <stdint.h>

#ifdef CONFIG_CORDIC_TRIGONOMETRY_FIXED
/* Use a local definition to avoid adding a dependency on <math.h> */
#define _M_PI		3.14159265358979323846	/* pi */
/**
 * \cordic_atan2_lookup_table = atan(2.^-(0:N-1)) N = 31/16
 * \CORDIC Gain is cordic_gain = prod(sqrt(1 + 2.^(-2*(0:31/16-1))))
 * \Inverse CORDIC Gain,inverse_cordic_gain = 1 / cordic_gain
 */
static const int32_t cordic_lookup[CORDIC_31B_TABLE_SIZE] = { 843314857, 497837829,
	263043837, 133525159, 67021687, 33543516, 16775851, 8388437, 4194283, 2097149,
	1048576, 524288, 262144, 131072, 65536, 32768, 16384, 8192, 4096, 2048, 1024,
	512, 256, 128, 64, 32, 16, 8, 4, 2, 1 };

/* 652032874 , deg = 69.586061*/
const int32_t cordic_sine_cos_lut_q29fl	 =  Q_CONVERT_FLOAT(1.214505869895220, 29);
/* 1686629713, deg = 90.000000	*/
const int32_t cordic_sine_cos_piovertwo_q30fl  = Q_CONVERT_FLOAT(_M_PI / 2, 30);
/* 421657428 , deg = 90.000000 */
const int32_t cord_sincos_piovertwo_q28fl  = Q_CONVERT_FLOAT(_M_PI / 2, 28);
/* 843314857,  deg = 90.000000	*/
const int32_t cord_sincos_piovertwo_q29fl  = Q_CONVERT_FLOAT(_M_PI / 2, 29);
/* arc trignometry constant*/
/**
 * \+----------+----------------------------------------+--------------------+-------------------+
 * \|values    |		Q_CONVERT_FLOAT		|(180/pi)*angleInRad |(pi/180)*angleInDeg|
 * \+----------+----------------------------------------+--------------------+-------------------+
 * \|379625062 | Q_CONVERT_FLOAT(1.4142135605216026, 28)|    81.0284683480568| 1.41421356052160  |
 * \|1073741824| Q_CONVERT_FLOAT(1.0000000000000000, 30)|    57.2957795130823| 1.00000000000000  |
 * \|843314856 | Q_CONVERT_FLOAT(1.5707963258028030, 29)|    89.9999999431572| 1.57079632580280  |
 * \|1686629713| Q_CONVERT_FLOAT(1.5707963267341256, 30)|    89.9999999965181| 1.57079632673413  |
 * \+----------+----------------------------------------+--------------------+-------------------+
 */
/* 379625062,  deg = 81.0284683480568475 or round(1.4142135605216026*2^28) */
const int32_t cord_arcsincos_q28fl  = Q_CONVERT_FLOAT(1.4142135605216026 / 2, 28);
/* 843314856 , deg = 89.9999999431572348 or round(pi/2*2^29) */
/* 843314857, Bug in calculation, so using hardcoded value */
const int32_t cord_arcsincos_piovertwo_q29fl  = Q_CONVERT_FLOAT(_M_PI / 2, 29);
/* 1073741824, deg = 57.2957795130823229 or round(1*2^30)*/
const int32_t cord_arcsincos_q30fl  = Q_CONVERT_FLOAT(1.0000000000000000, 30);

/**
 * \CORDIC-based approximation of sine, cosine and complex exponential
 */
void cordic_approx(int32_t th_rad_fxp, int32_t a_idx, int32_t *sign, int32_t *b_yn, int32_t *xn,
		   int32_t *th_cdc_fxp)
{
	int32_t b_idx;
	int32_t xtmp;
	int32_t ytmp;
	*sign = 1;
	/* Addition or subtraction by a multiple of pi/2 is done in the data type
	 * of the input. When the fraction length is 29, then the quantization error
	 * introduced by the addition or subtraction of pi/2 is done with 29 bits of
	 * precision.Input range of cordicsin must be in the range [-2*pi, 2*pi),
	 * a signed type with fractionLength = wordLength-4 will fit this range
	 * without overflow.Increase of fractionLength makes the addition or
	 * subtraction of a multiple of pi/2 more precise
	 */
	if (th_rad_fxp > cord_sincos_piovertwo_q28fl) {
		if ((th_rad_fxp - cord_sincos_piovertwo_q29fl) <= cord_sincos_piovertwo_q28fl) {
			th_rad_fxp -= cord_sincos_piovertwo_q29fl;
			*sign  = -1;
		} else {
			th_rad_fxp -= cordic_sine_cos_piovertwo_q30fl;
		}
	} else if (th_rad_fxp < -cord_sincos_piovertwo_q28fl) {
		if ((th_rad_fxp + cord_sincos_piovertwo_q29fl) >= -cord_sincos_piovertwo_q28fl) {
			th_rad_fxp += cord_sincos_piovertwo_q29fl;
			*sign  = -1;
		} else {
			th_rad_fxp += cordic_sine_cos_piovertwo_q30fl;
		}
	}

	th_rad_fxp <<= 2;
	*b_yn = 0;
	*xn = cordic_sine_cos_lut_q29fl;
	xtmp = cordic_sine_cos_lut_q29fl;
	ytmp = 0;

	/* Calculate the correct coefficient values from rotation angle.
	 * Find difference between the coefficients from the lookup table
	 * and those from the calculation
	 */
	for (b_idx = 0; b_idx < a_idx; b_idx++) {
		if (th_rad_fxp < 0) {
			th_rad_fxp += cordic_lookup[b_idx];
			*xn += ytmp;
			*b_yn -= xtmp;
		} else {
			th_rad_fxp -= cordic_lookup[b_idx];
			*xn -= ytmp;
			*b_yn += xtmp;
		}
		xtmp = *xn >> (b_idx + 1);
		ytmp = *b_yn >> (b_idx + 1);
	}
	/* Q2.30 format -sine,cosine*/
	*th_cdc_fxp = th_rad_fxp;
}

/**
 * \CORDIC-based approximation for inverse sine, cosine
 * Arguments	: int32_t realvalue
 *		  int16_t numiters
 *		  const int32_t LUT[CORDIC_30B_ITABLE_SIZE][30]
 * Return Type	: int32_t
 */
int32_t is_scalar_cordic(int32_t realvalue, int16_t numiters,
			 const int32_t LUT[CORDIC_30B_ITABLE_SIZE], cordic_cfg type)
{
	int32_t b_sincosvalue;
	int32_t b_i;
	int32_t i;
	int32_t x = 0;
	int32_t xdshift;
	int32_t xshift;
	int32_t y = 0;
	int32_t ydshift;
	int32_t yshift;
	int32_t z = 0;
	int16_t j;
	int16_t k;

	b_sincosvalue = (realvalue >> 1U);
	if (type == EN_32B_CORDIC_ASINE) {
		/* Initialize the variables for the cordic iteration
		 * angles less than pi/4, we initialize (x,y) along the x-axis.
		 * angles greater than or equal to pi/4, we initialize (x,y)
		 * along the y-axis. This improves the accuracy of the algorithm
		 * near the edge of the domain of convergence
		 */
		if (b_sincosvalue > cord_arcsincos_q28fl) {
			x = 0;
			y = cord_arcsincos_q30fl;
			z = 843314856;
		} else {
			x = cord_arcsincos_q30fl;
			y = 0;
			z = 0;
		}
	}
	if (type == EN_32B_CORDIC_ACOSINE) {
		/* Initialize the variables for the cordic iteration
		 * angles less than pi/4, we initialize (x,y) along the x-axis.
		 * angles greater than or equal to pi/4, we initialize (x,y)
		 * along the y-axis. This improves the accuracy of the algorithm
		 * near the edge of the domain of convergence
		 */
		if (b_sincosvalue < cord_arcsincos_q28fl) {
			x = 0;
			y = cord_arcsincos_q30fl;
			z = 843314856;
		} else {
			x = cord_arcsincos_q30fl;
			y = 0;
			z = 0;
		}
	}
	b_sincosvalue <<= 1U;
	/* DCORDIC(Double CORDIC) algorithm */
	/* Double iterations method consists in the fact that unlike the classical */
	/* CORDIC method,where the iteration step value changes EVERY time, i.e. on */
	/* each iteration, in the double iteration method, the iteration step value */
	/* is repeated twice and changes only through one iteration */
	i = numiters - 1;
	for (b_i = 0; b_i < i; b_i++) {
		j = (int16_t)((b_i + 1) << 1U);
		if (j >= 31)
			j = 31;

		if (b_i < 31)
			k = (int16_t)b_i;
		else
			k = 31;

		xshift = (x >>	(k));
		xdshift = (x >> (j));
		yshift = (y >>	(k));
		ydshift = (y >> (j));

		if (type == EN_32B_CORDIC_ASINE) {
			/* Do nothing if x currently equals the target value. Allowed for
			 * double rotations algorithms, as it is equivalent to rotating by
			 * the same angle in opposite directions sequentially. Accounts for
			 * the scaling effect of CORDIC Pseudo-rotations as well.
			 */
			if (y == b_sincosvalue) {
				x += xdshift;
				y += ydshift;
			} else if (((y >= b_sincosvalue) && (x >= 0)) ||
				  ((y < b_sincosvalue) && (x < 0))) {
				x = (x - xdshift) + yshift;
				y = (y - ydshift) - xshift;
				z -= LUT[b_i];
			} else {
				x = (x - xdshift) - yshift;
				y = (y - ydshift) + xshift;
				z += LUT[b_i];
			}
		}
		if (type == EN_32B_CORDIC_ACOSINE) {
			/* Do nothing if x currently equals the target value. Allowed for
			 * double rotations algorithms, as it is equivalent to rotating by
			 * the same angle in opposite directions sequentially. Accounts for
			 * the scaling effect of CORDIC Pseudo-rotations as well.
			 */
			if (x == b_sincosvalue) {
				x += xdshift;
				y += ydshift;
			} else if (((x > b_sincosvalue) && (y >= 0)) ||
				  ((x < b_sincosvalue) && (y < 0))) {
				x = (x - xdshift) - yshift;
				y = (y - ydshift) + xshift;
				z += LUT[b_i];
			} else {
				x = (x - xdshift) + yshift;
				y = (y - ydshift) - xshift;
				z -= LUT[b_i];
			}
		}

		b_sincosvalue += (b_sincosvalue >> ((uint32_t)j));
	}
	if (z < 0)
		z = -z;

	return z;
}

/**
 * approximated complex result
 * Arguments	: int32_t sign
 *		  int32_t b_yn
 *		  int32_t xn
 *		  enum    type
 * Return Type	: struct
 */
struct cmpx cmpx_cexp(int32_t sign, int32_t b_yn, int32_t xn, cordic_cfg type)
{
	struct cmpx cexp;

	cexp.re = sign * xn;
	cexp.im = sign * b_yn;
	/*convert Q2.30 to Q17.15*/
	if (type == EN_16B_CORDIC_CEXP) {
		cexp.re =  sat_int16(Q_SHIFT_RND((cexp.re), 30, 15));
		cexp.im = sat_int16(Q_SHIFT_RND((cexp.im), 30, 15));
	}
	return cexp;
}

#endif
