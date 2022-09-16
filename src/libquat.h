/*
 * CDDL HEADER START
 *
 * The contents of this file are subject to the terms of the
 * Common Development and Distribution License, Version 1.0 only
 * (the "License").  You may not use this file except in compliance
 * with the License.
 *
 * You can obtain a copy of the license in the file COPYING
 * or http://www.opensource.org/licenses/CDDL-1.0.
 * See the License for the specific language governing permissions
 * and limitations under the License.
 *
 * When distributing Covered Code, include this CDDL HEADER in each
 * file and include the License file COPYING.
 * If applicable, add the following below this CDDL HEADER, with the
 * fields enclosed by brackets "[]" replaced with your own identifying
 * information: Portions Copyright [yyyy] [name of copyright owner]
 *
 * CDDL HEADER END
 */
/*
 * Copyright 2020 Saso Kiselkov. All rights reserved.
 */

#ifndef	_LIBQUAT_H_
#define	_LIBQUAT_H_

#include <acfutils/geom.h>
#include <mathc.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * This library makes it easier to use the mathc library quaternion
 * functions for keeping track of rotations. Internally, libquat uses
 * the (X, Y, Z, W) quaternion format, which matches the mathc library.
 * However, X-Plane uses the (W, X, Y, Z) format. To convert between
 * the two formats, use the "quat_from_xp" and "quat_to_xp" functions.
 *
 * Also note that axis assignments in OpenGL space and quaternion space
 * are different. Use QUAT_AXIS_*_GL to pick out the correct coordinate
 * to convert into OpenGL vectors. In either case, however, angles are
 * in right-hand rule (i.e. with the right thumb along the axis, the
 * fingers are pointing in the direction of positive rotation angles).
 */

#if defined(MATHC_USE_UNIONS)
#define	QUAT_VEC3(_x, _y, _z) ((struct vec3){.v = { [0] = (mfloat_t)(_x), [1] = (mfloat_t)(_y), [2] = (mfloat_t)(_z)}})
#define QUAT_VEC4(_x, _y, _z, _w) ((struct quat){.v = { [0] = (mfloat_t)(_x), [1] = (mfloat_t)(_y), [2] = (mfloat_t)(_z), [3] = (mfloat_t)(_w)}})
#else
#define QUAT_VEC3(_x, _y, _z) ((struct vec3){.x = (mfloat_t)(_x), .y = (mfloat_t)(_y), .z = (mfloat_t)(_z)})
#define QUAT_VEC4(_x, _y, _z, _w) ((struct quat){.x = (mfloat_t)(_x), .y = (mfloat_t)(_y), .z = (mfloat_t)(_z), .w = (mfloat_t)(_w)})
#endif

/*
 * Unit vectors along the OpenGL X, Y and Z axes in libquat internal format.
 */
#define	QUAT_AXIS_X_GL		QUAT_VEC3(0, 1, 0)
#define	QUAT_AXIS_Y_GL		QUAT_VEC3(0, 0, -1)
#define	QUAT_AXIS_Z_GL		QUAT_VEC3(-1, 0, 0)

#define	IS_NULL_QUAT(q)		(!isfinite((q).w))
#define	NULL_QUAT		((struct quat){{{NAN, NAN, NAN, NAN}}})

struct quat quat_hamil_prod(struct quat p, struct quat q);
struct quat quat_local2ecmigl(geo_pos2_t refpt, mfloat_t ref_time);
struct quat quat_ecmigl2local(geo_pos2_t refpt, mfloat_t ref_time);
struct quat quat_rot_rel(struct quat q1, struct quat q2);
struct quat quat_rot_concat(struct quat from, struct quat delta);
void quat_to_axis_angle(struct quat q, struct vec3 *axis, mfloat_t *angle);
void quat_to_axis_gl_angle(struct quat q, struct vec3 *axis, mfloat_t *angle);
struct quat quat_from_euler(mfloat_t psi, mfloat_t theta, mfloat_t phi)
    PURE_ATTR;
struct quat quat_from_euler_deg(mfloat_t psi, mfloat_t theta, mfloat_t phi)
    PURE_ATTR;
void quat_to_euler(struct quat q, mfloat_t *psi, mfloat_t *theta,
    mfloat_t *phi);
void quat_to_euler_deg(struct quat q, mfloat_t *psi, mfloat_t *theta,
    mfloat_t *phi);
struct quat quat_from_xp(struct quat xp_q) PURE_ATTR;
struct quat quat_to_xp(struct quat q) PURE_ATTR;

/*
 * Converts a vect3l_t in OpenGL space to a quaternion suitable for
 * transformation using Hamilton products.
 */
static inline struct quat
quat_from_vect3l_gl(vect3l_t v)
{
	struct quat q = {};

	q.x = -v.z;
	q.y = v.x;
	q.z = -v.y;

	return (q);
}

/*
 * Converts a quat back to a vect3l_t in OpenGL space.
 */
static inline vect3l_t
quat_to_vect3l_gl(struct quat q)
{
	return (VECT3L(q.y, -q.z, -q.x));
}

/*
 * Converts a quat to a vect3l_t in North-East-Down orientation.
 * In reality it just returns the XYZ components, since quaternions
 * are already in NED notation.
 */
static inline vect3l_t
quat_to_vect3l_NED(struct quat q)
{
	return (VECT3L(q.x, q.y, q.z));
}

/*
 * Shorthand for quat_inverse that allows using the assignment notation
 *	struct quat q = quat_inv(p);
 * Instead of having to declare the quaternion ahead of time and then
 * passing both quaternion vector fields by reference.
 */
static inline struct quat
quat_inv(struct quat q)
{
	struct quat p;
#if defined(MATHC_USE_UNIONS)
	quat_inverse(p.v, q.v);
#else
    quat_inverse((mfloat_t*) &p, (mfloat_t*) &q);
#endif
	return (p);
}

void quat_print(struct quat q);

void
quat_print_named(char *name, struct quat q);

void 
_euler(struct quat q);
void quat_print_axis_angle(struct quat q);

/*
 * North-East-Down coordinates
 */
typedef struct {
	mfloat_t N;
	mfloat_t E;
	mfloat_t D;
} NED_t;

static inline struct quat
NED2quat(NED_t ned)
{
#if defined(MATHC_USE_UNIONS)
	return ((struct quat){{{ ned.N, ned.E, ned.D }}});
#else
    return ((struct quat) {.x = ned.N, .y = ned.E, .z = ned.D, .w = 0});
#endif
}

static inline NED_t
quat2NED(struct quat q)
{
	return ((NED_t){ q.x, q.y, q.z });
}

#define	NULL_NED		((NED_t){NAN, NAN, NAN})
#define	IS_NULL_NED(NED)	(isnan((NED).N))

#ifdef __cplusplus
}
#endif

#endif	/* _LIBQUAT_H_ */
