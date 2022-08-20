#![allow(dead_code)]

use core::{
    f64::consts::TAU,
    ops::{Add, Mul, Sub},
};

const EPS: f64 = 0.0000001;

#[derive(Clone, Copy, Default, Debug)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Mul<f64> for Vec3 {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl Add<Self> for Vec3 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl Sub<Self> for Vec3 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl Vec3 {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn new_zero() -> Self {
        Self {
            x: 0.,
            y: 0.,
            z: 0.,
        }
    }

    /// Returns the vector magnitude.
    pub fn magnitude(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }

    /// Returns the normalised version of the vector
    pub fn to_normalized(self) -> Self {
        let mag_recip = 1. / self.magnitude();
        self * mag_recip
    }

    /// Calculate the dot product
    pub fn dot(&self, rhs: Self) -> f64 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    /// Calculate the cross product.
    pub fn cross(&self, rhs: Self) -> Self {
        Self {
            x: self.y * rhs.z - self.z * rhs.y,
            y: self.z * rhs.x - self.x * rhs.z,
            z: self.x * rhs.y - self.y * rhs.x,
        }
    }

    /// Project a vector onto a plane defined by its normal vector. Assumes self and `plane_norm`
    /// are unit vectors.
    pub fn project_to_plane(&self, plane_norm: Self) -> Self {
        *self - plane_norm * self.dot(plane_norm)
    }
}

#[derive(Clone, Copy, Default, Debug)]
pub struct Quaternion {
    pub w: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Mul<Self> for Quaternion {
    type Output = Self;

    /// Multiply a quaternion by another quaternion. This can be used to... todo
    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            w: self.w * rhs.w - self.x * rhs.x - self.y * rhs.y - self.z * rhs.z,
            x: self.w * rhs.x + self.x * rhs.w + self.y * rhs.z - self.z * rhs.y,
            y: self.w * rhs.y - self.x * rhs.z + self.y * rhs.w + self.z * rhs.x,
            z: self.w * rhs.z + self.x * rhs.y - self.y * rhs.x + self.z * rhs.w,
        }
    }
}

impl Mul<Vec3> for Quaternion {
    type Output = Self;

    /// Returns the multiplication of a Quaternion with a vector.  This is a
    /// normal Quaternion multiplication where the vector is treated as a
    /// Quaternion with a W element value of zero. The Quaternion is post-
    /// multiplied by the vector.
    fn mul(self, rhs: Vec3) -> Self::Output {
        Self {
            w: -self.x * rhs.x - self.y * rhs.y - self.z * rhs.z,
            x: self.w * rhs.x + self.y * rhs.z - self.z * rhs.y,
            y: self.w * rhs.y - self.x * rhs.z + self.z * rhs.x,
            z: self.w * rhs.z + self.x * rhs.y - self.y * rhs.x,
        }
    }
}

impl Mul<f64> for Quaternion {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            w: self.w * rhs,
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl Quaternion {
    pub fn new_identity() -> Self {
        Self {
            w: 1.,
            x: 0.,
            y: 0.,
            z: 0.,
        }
    }

    /// Create the quaternion that creates the shortest (great circle) rotation from vec0
    /// to vec1.
    pub fn from_unit_vecs(v0: Vec3, v1: Vec3) -> Self {
        const ONE_MINUS_EPS: f64 = 1.0 - 2.0 * EPS;

        let dot = v0.dot(v1);
        if dot > ONE_MINUS_EPS {
            return Self::new_identity();
        } else if dot < -ONE_MINUS_EPS {
            // Rotate along any orthonormal vec to vec1 or vec2 as the axis.
            return Self::from_axis_angle(Vec3::new(1., 0., 0.).cross(v0), TAU / 2.);
        }

        let w = 1. + dot;
        let v = v0.cross(v1);

        (Self {
            w,
            x: v.x,
            y: v.y,
            z: v.z,
        })
        .to_normalized()
    }

    pub fn from_euler(phi: f64, psi: f64, theta: f64) -> Self {
        let cy = (theta * 0.5).cos();
        let sy = (theta * 0.5).sin();
        let cp = (psi * 0.5).cos();
        let sp = (psi * 0.5).sin();
        let cr = (phi * 0.5).cos();
        let sr = (phi * 0.5).sin();

        Self {
            w: cr * cp * cy + sr * sp * sy,
            x: sr * cp * cy - cr * sp * sy,
            y: cr * sp * cy + sr * cp * sy,
            z: cr * cp * sy - sr * sp * cy,
        }
    }

    /// Convert this quaternion to Euler angles.
    /// https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
    pub fn to_euler(&self) -> (f64, f64, f64) {
        // roll, pitch, yaw, (x, y, z axes)
        // roll (x-axis rotation)
        let sinr_cosp = 2. * (self.w * self.x + self.y * self.z);
        let cosr_cosp = 1. - 2. * (self.x * self.x + self.y * self.y);

        let roll = sinr_cosp.atan2(cosr_cosp);

        // pitch (y-axis rotation)
        let sinp = 2. * (self.w * self.y - self.z * self.x);
        let pitch = if sinp.abs() >= 1. {
            (TAU / 4.).copysign(sinp) // use 90 degrees if out of range
        } else {
            sinp.asin()
        };

        // yaw (z-axis rotation)
        let siny_cosp = 2. * (self.w * self.z + self.x * self.y);
        let cosy_cosp = 1. - 2. * (self.y * self.y + self.z * self.z);
        let yaw = siny_cosp.atan2(cosy_cosp);

        (roll, pitch, yaw)
    }

    // /// Creates an orientation that point towards a vector, with a given up direction defined.
    // pub fn from_vec_direction(dir: Vec3, up: Vec3) -> Self {
    //     let forward_vector = dir;
    //
    //     let forward = Vec3::new(0., 0., 1.);
    //
    //     let dot = forward.dot(forward_vector);
    //
    //     if (dot - (-1.0)).abs() < 0.000001 {
    //         // return Self: { x:  Quaternion(Vector3.up.x, Vector3.up.y, Vector3.up.z, 3.1415926535897932f);
    //         Self::new_identity(); // todo! adapt the above.
    //
    //     }
    //     if (dot - (1.0)).abs() < 0.000001 {
    //         return Self::new_identity();
    //     }
    //
    //     let rot_angle = dot.acos();
    //     let rot_axis = forward.cross(forward_vector).to_normalized();
    //
    //     Self::from_axis_angle(rot_axis, rot_angle)
    // }

    pub fn inverse(self) -> Self {
        Self {
            w: self.w,
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }

    /// Rotate a vector using this quaternion. Note that our multiplication Q * v
    /// operation is effectively quaternion multiplication, with a quaternion
    /// created by a vec with w=0.
    pub fn rotate_vec(self, vec: Vec3) -> Vec3 {
        (self * vec * self.inverse()).to_vec()
    }

    /// Create a rotation quaternion from an axis and angle.
    /// `axis` mus be normalized.
    pub fn from_axis_angle(axis: Vec3, angle: f64) -> Self {
        // Here we calculate the sin( theta / 2) once for optimization
        let factor = (angle / 2.).sin();

        Self {
            // Calcualte the w value by cos( theta / 2 )
            w: (angle / 2.).cos(),
            // Calculate the x, y and z of the quaternion
            x: axis.x * factor,
            y: axis.y * factor,
            z: axis.z * factor,
        }
    }

    /// Convert to a 3D vector, discarding `w`.
    pub fn to_vec(self) -> Vec3 {
        Vec3 {
            x: self.x,
            y: self.y,
            z: self.z,
        }
    }

    /// Returns the vector magnitude.
    pub fn magnitude(&self) -> f64 {
        (self.w.powi(2) + self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }

    /// Returns the normalised version of the vector
    pub fn to_normalized(self) -> Self {
        let mag_recip = 1. / self.magnitude();
        self * mag_recip
    }
}

/// Calculate the determinate of a matrix defined by its columns.
/// We use this for determining the full 0 - tau angle between bonds.
#[rustfmt::skip]
pub fn det_from_cols(c0: Vec3, c1: Vec3, c2: Vec3) -> f64 {
    c0.x * c1.y * c2.z +
    c1.x * c2.y * c0.z +
    c2.x * c0.y * c1.z -
    c0.x * c2.y * c1.z -
    c1.x * c0.y * c2.z -
    c2.x * c1.y * c0.z
}
