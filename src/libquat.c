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

#include <stdio.h>

#include "libquat.h"

#include <acfutils/assert.h>
#include <acfutils/perf.h>
#include <acfutils/log.h>


CTASSERT(sizeof (mfloat_t) == sizeof (long double));

/*
 * Calculates the Hamilton product of quaternions `p' and `q':
 * H(p,q) = q.p.q^-1
 */
struct quat
quat_hamil_prod(struct quat p, struct quat q)
{
	struct quat q_inv, out;

	ASSERT(!IS_NULL_QUAT(p));
	ASSERT(!IS_NULL_QUAT(q));

#if defined(MATHC_USE_UNIONS)
	quat_inverse(q_inv.v, q.v);
	quat_multiply(out.v, q.v, p.v);
	quat_multiply(out.v, out.v, q_inv.v);
#else
	quat_inverse(&q_inv, &q);
	quat_multiply(&out, &q, &p);
	quat_multiply(&out, &out, &q_inv);
#endif

	return (out);
}

/*
 * Constructs a quaternion to translate from X-Plane OpenGL local coordinates
 * with origin at `refpt' to the Earth-Centered-modified-Inertial frame in
 * OpenGL coordinates at `ref_time'.
 *
 * X-Plane OpenGL local coordinates are a euclidian coordinate space with
 * three axes at right angles to each other. The coordinate system origin is
 * at a given reference latitude/longitude point, at elevation 0 meters
 * above the Earth ellipsoid. The X axis is parallel with the Earth's surface
 * at this reference point and increases towards to the east. The Y axis is
 * perpendicular to the Earth's surface at the reference point and increases
 * upward. The Z axis is parallel to the Earth's surface and increases to the
 * south. All coordinates are in meters from the origin point.
 *
 * ECmIGL coordinates work as follows. ECmI is a euclidian coordinate space
 * composed of three axes (X, Y & Z) at right angles to each other. The
 * origin of the coordinate system is centered on the Earth's geometric
 * center of rotation. Being an inertial frame, the axes only line up with
 * certain reference points on the Earth's surface at ref_time=0. Thereafter,
 * the Earth rotates counter-clockwise in the coordinate system, which the
 * axes of the coordinate system remain aligned with the fixed background
 * stars. `ref_time' defines the rotation offset of the Earth in seconds
 * from the starting orientation.
 * At `ref_time=0', the X axis goes from the origin through the latitude=0
 * and longitude=90 degree point on the Earth. The Y axis goes
 * from the origin to the latitude=90 (north pole) point. The Z axis goes
 * from the origin to the latitude=0 and longitude=0 point. All coordinates
 * are in meters from the origin point.
 */
struct quat
quat_local2ecmigl(geo_pos2_t refpt, mfloat_t ref_time)
{
	struct quat lat_q, lon_q, out_q;

	ASSERT(isfinite(ref_time));
	ASSERT(!IS_NULL_GEO_POS(refpt));
#if defined(MATHC_USE_UNIONS)
	
	struct vec3 open_gl_x_axis_v3 = {0.0, 1.0, 0.0};
	struct vec3 open_gl_y_axis_v3 = {0.0, 0.0, -1.0};

	quat_from_axis_angle(lat_q.v, open_gl_x_axis_v3.v, DEG2RAD(90 - refpt.lat));
	quat_from_axis_angle(lon_q.v, open_gl_y_axis_v3.v, DEG2RAD(refpt.lon + EARTH_ROT_RATE * ref_time));

	// logMsg("[DEBUG] lon_q is %lf, %lf, %lf, %lf", (double) lon_q.x, (double) lon_q.y, 
 //                                                                (double) lon_q.z, (double) lon_q.w);

 //    logMsg("[DEBUG] lat_q is %lf, %lf, %lf, %lf", (double) lat_q.x, (double) lat_q.y, 
 //                                                                (double) lat_q.z, (double) lat_q.w);

	quat_multiply(out_q.v, lon_q.v, lat_q.v);

    // logMsg("[DEBUG] out_q is %lf, %lf, %lf, %lf", (double) out_q.x, (double) out_q.y, 
    //                                                            (double) out_q.z, (double) out_q.w);

#else 
	out_q = (struct quat) {.x = 0, .y = 0, .z = 0, .w = 0};
	quat_from_axis_angle((mfloat_t*) &lat_q, (mfloat_t*) &QUAT_AXIS_X_GL,
	    DEG2RAD(90 - refpt.lat));
	quat_from_axis_angle( (mfloat_t*)  &lon_q,(mfloat_t*) &QUAT_AXIS_Y_GL,
	    DEG2RAD(refpt.lon + EARTH_ROT_RATE * ref_time));

	// logMsg("[DEBUG] lon_q is %lf, %lf, %lf, %lf", (double) lon_q.x, (double) lon_q.y, 
 //                                                                (double) lon_q.z, (double) lon_q.w);

 //    logMsg("[DEBUG] lat_q is %lf, %lf, %lf, %lf", (double) lat_q.x, (double) lat_q.y, 
 //                                                                (double) lat_q.z, (double) lat_q.w);

 //    logMsg("[DEBUG] out_q is %lf, %lf, %lf, %lf", (double) out_q.x, (double) out_q.y, 
 //                                                                (double) out_q.z, (double) out_q.w);

	quat_multiply((mfloat_t*) &out_q, (mfloat_t*) &lon_q, (mfloat_t*) &lat_q);

	// logMsg("[DEBUG] out_q is %lf, %lf, %lf, %lf", (double) out_q.x, (double) out_q.y, 
    //                                                            (double) out_q.z, (double) out_q.w);
#endif

	return (out_q);
}

/*
 * Constructs inverse quaternion to translate from ECmI to X-Plane OpenGL
 * local coordinates. See quat_local2ecmigl for a description of the
 * respective coordinate systems.
 */
struct quat
quat_ecmigl2local(geo_pos2_t refpt, mfloat_t ref_time)
{
	struct quat q, q_inv;

	ASSERT(isfinite(ref_time));
	ASSERT(!IS_NULL_GEO_POS(refpt));

	q = quat_local2ecmigl(refpt, ref_time);

#if defined(MATHC_USE_UNIONS)
	quat_inverse(q_inv.v, q.v);
#else 
	quat_inverse((mfloat_t*) &q_inv, (mfloat_t*) &q);
#endif

	return (q_inv);
}

/*
 * Constructs a relative rotation quaternion that expresses the relative
 * rotation from `q1' to `q2'. This can subsequently be concatenated onto
 * `q1' to obtain `q2'.
 */
struct quat
quat_rot_rel(struct quat q1, struct quat q2)
{
	struct quat q1_inv, d_quat;

	ASSERT(!IS_NULL_QUAT(q1));
	ASSERT(!IS_NULL_QUAT(q2));

#if defined(MATHC_USE_UNIONS)
	quat_inverse(q1_inv.v, q1.v);
	quat_multiply(d_quat.v, q2.v, q1_inv.v);
#else
	quat_inverse((mfloat_t*) &q1_inv, (mfloat_t*) &q1);
	quat_multiply((mfloat_t*) &d_quat, (mfloat_t*) &q2, (mfloat_t*) &q1_inv);
#endif

	return (d_quat);
}

/*
 * Concatenates a relative rotation quaternion `delta' onto `from'.
 * This is equivalent to quaternion multiplication of `delta' x `from'.
 */
struct quat
quat_rot_concat(struct quat from, struct quat delta)
{
	struct quat to;

	ASSERT(!IS_NULL_QUAT(from));
	ASSERT(!IS_NULL_QUAT(delta));

#if defined(MATHC_USE_UNIONS)
	quat_multiply(to.v, delta.v, from.v);
#else 
	quat_multiply((mfloat_t*) &to, (mfloat_t*) &delta, (mfloat_t*) &from);
#endif

	return (to);
}

/*
 * Converts a quaternion into an axis + angle (in radians) representation.
 * For any argument you don't wish to receive, pass NULL.
 */
void
quat_to_axis_angle(struct quat q, struct vec3 *axis, mfloat_t *angle)
{
	mfloat_t half, s;

	ASSERT(!IS_NULL_QUAT(q));

	half = acosl(q.w);
	s = sinl(half);
	/*
	 * Zero roration quaternion has an arbitrary axis, so just pick one.
	 */
	if (s == 0 || isnan(s)) {
		if (axis != NULL)
			*axis = QUAT_AXIS_X_GL;
		if (angle != NULL)
			*angle = 0;
		return;
	}
	if (axis != NULL) {
		axis->x = q.x / s;
		axis->y = q.y / s;
		axis->z = q.z / s;
	}
	if (angle != NULL)
		*angle = 2 * half;
}

/*
 * Same as quat_to_axis_angle, but converts the axis into an axis in
 * OpenGL coordinates.
 */
void
quat_to_axis_gl_angle(struct quat q, struct vec3 *axis, mfloat_t *angle)
{
	ASSERT(!IS_NULL_QUAT(q));
	quat_to_axis_angle(q, axis, angle);
	if (axis != NULL)
		*axis = QUAT_VEC3(axis->y, -axis->z, -axis->x);
}

/*
 * Creates a quaternion from a set of Euler angles (yaw, pitch, roll)
 * in radians.
 */
struct quat
quat_from_euler(mfloat_t psi, mfloat_t theta, mfloat_t phi)
{
	struct quat q;
	mfloat_t sin_psi = sinl(psi / 2), cos_psi = cosl(psi / 2);
	mfloat_t sin_theta = sinl(theta / 2), cos_theta = cosl(theta / 2);
	mfloat_t sin_phi = sinl(phi / 2), cos_phi = cosl(phi / 2);

	ASSERT(isfinite(psi));
	ASSERT(isfinite(theta));
	ASSERT(isfinite(phi));
	q.w =  cos_psi * cos_theta * cos_phi + sin_psi * sin_theta * sin_phi;
	q.x =  cos_psi * cos_theta * sin_phi - sin_psi * sin_theta * cos_phi;
	q.y =  cos_psi * sin_theta * cos_phi + sin_psi * cos_theta * sin_phi;
	q.z = -cos_psi * sin_theta * sin_phi + sin_psi * cos_theta * cos_phi;

	return (q);
}

/*
 * Creates a quaternion from a set of Euler angles (yaw, pitch, roll)
 * in degrees.
 */
struct quat
quat_from_euler_deg(mfloat_t psi, mfloat_t theta, mfloat_t phi)
{
	return (quat_from_euler(DEG2RAD(psi), DEG2RAD(theta), DEG2RAD(phi)));
}

/*
 * Converts a quaternion into Euler angles (yaw, pitch, roll) in radians.
 * The pointer output arguments are optional - any that are left NULL are
 * not populated.
 */
void
quat_to_euler(struct quat q, mfloat_t *psi, mfloat_t *theta, mfloat_t *phi)
{
	mfloat_t sq_w = q.w * q.w;
	mfloat_t sq_x = q.x * q.x;
	mfloat_t sq_y = q.y * q.y;
	mfloat_t sq_z = q.z * q.z;

	ASSERT(!IS_NULL_QUAT(q));
	if (psi != NULL) {
		*psi = atan2l(2.0 * (q.x * q.y + q.z * q.w),
		    (sq_x - sq_y - sq_z + sq_w));
	}
	if (theta != NULL)
		*theta = asinl(-2.0 * (q.x * q.z - q.y * q.w));
	if (phi != NULL) {
		*phi = atan2l(2.0 * (q.y * q.z + q.x * q.w),
		    (-sq_x - sq_y + sq_z + sq_w));
	}
}

/*
 * Converts a quaternion into Euler angles (yaw, pitch, roll) in degrees.
 * The pointer output arguments are optional - any that are left NULL are
 * not populated.
 */
void
quat_to_euler_deg(struct quat q, mfloat_t *psi, mfloat_t *theta, mfloat_t *phi)
{
	quat_to_euler(q, psi, theta, phi);
	if (psi != NULL)
		*psi = RAD2DEG(*psi);
	if (theta != NULL)
		*theta = RAD2DEG(*theta);
	if (phi != NULL)
		*phi = RAD2DEG(*phi);
}

/*
 * Transforms a quaternion in X-Plane format (W, X, Y, Z) into libquat
 * internal format (X, Y, Z, W).
 */
struct quat
quat_from_xp(struct quat xp_q)
{
#if defined(MATHC_USE_UNIONS)
	struct quat q = { .v = { xp_q.v[1], xp_q.v[2], xp_q.v[3], xp_q.v[0] } };
#else
	struct quat q = { .x = xp_q.y, .y = xp_q.z, .z = xp_q.w, .w = xp_q.x };
#endif

	return (q);
}

/*
 * Transforms a quaternion in libquat internal format (X, Y, Z, W) into
 * the format that X-Plane uses in its rotation quaternion (W, X, Y, Z).
 */
struct quat
quat_to_xp(struct quat q)
{

#if defined(MATHC_USE_UNIONS)
	struct quat xp_q = { .v = { q.v[3], q.v[0], q.v[1], q.v[2] } };
#else 
	struct quat xp_q = { .x = q.w, .y = q.x, .z = q.y, .w = q.z };
#endif

	return (xp_q);
}

void
quat_print(struct quat q)
{
#ifndef	__MINGW32__
	logMsg(".w=%11Lf .x=%11Lf .y=%11Lf .z=%11Lf\n", q.w, q.x, q.y, q.z);
#else	/* defined(__MINGW32__) */
	logMsg(".w=%11f .x=%11f .y=%11f .z=%11f\n",
	    (double)q.w, (double)q.x, (double)q.y, (double)q.z);
#endif	/* defined(__MINGW32__) */
}

void
quat_print_named(char *name, struct quat q)
{
#ifndef	__MINGW32__
	logMsg("Quat %s : .w=%11Lf .x=%11Lf .y=%11Lf .z=%11Lf\n", name, q.w, q.x, q.y, q.z);
#else	/* defined(__MINGW32__) */
	logMsg("Quat %s : .w=%11f .x=%11f .y=%11f .z=%11f\n",
	    name, (double)q.w, (double)q.x, (double)q.y, (double)q.z);
#endif	/* defined(__MINGW32__) */
}


/*
 * Debug prints a quaternion transformed into Euler angle space
 * (yaw, pitch, roll), in degrees.
 */
void
quat_print_euler(struct quat q)
{
	mfloat_t psi, theta, phi;

	ASSERT(!IS_NULL_QUAT(q));
	quat_to_euler(q, &psi, &theta, &phi);
#ifndef	__MINGW32__
	logMsg("psi: %Lf  theta: %Lf  phi: %Lf",
	    RAD2DEG(psi), RAD2DEG(theta), RAD2DEG(phi));
#else	/* defined(__MINGW32__) */
	logMsg("psi: %f  theta: %f  phi: %f",
	    (double)RAD2DEG(psi), (double)RAD2DEG(theta), (double)RAD2DEG(phi));
#endif	/* defined(__MINGW32__) */
}

/*
 * Debug prints a quaternion as an axis + angle. The axis is automatically
 * translated into OpenGL coordinates.
 */
void
quat_print_axis_angle(struct quat q)
{
	struct vec3 axis;
	mfloat_t angle;

	ASSERT(!IS_NULL_QUAT(q));
	quat_to_axis_angle(q, &axis, &angle);
#ifndef	__MINGW32__
	logMsg("%11Lf %11Lf %11Lf %11Lf deg", axis.y, -axis.z, -axis.x,
	    RAD2DEG(angle));
#else	/* defined(__MINGW32__) */
	logMsg("%11f %11f %11f %11f deg", (double)axis.y, (double)-axis.z,
	    (double)-axis.x, (double)RAD2DEG(angle));
#endif	/* defined(__MINGW32__) */
}