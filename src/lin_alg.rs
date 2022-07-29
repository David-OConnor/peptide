use core::ops::{Add, Mul};

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

impl Vec3 {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
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
    /// normal Quaternion multiplication where the vector is treated a
    /// Quaternion with a W element value of zero.  The Quaternion is post-
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

impl Quaternion {
    pub fn new_identity() -> Self {
        Self {
            w: 1.,
            x: 0.,
            y: 0.,
            z: 0.,
        }
    }

    pub fn inverse(self) -> Self {
        Self {
            w: self.w,
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }

    /// Rotate a vector using this quaternion.
    pub fn rotate_vec(self, vec: Vec3) -> Vec3 {
        (self * vec * self.inverse()).to_vec()
    }

    /// Create a rotation quaternion from an axis and angle.
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
}
//
// /// A point in cartesian coordinates
// #[derive(Clone, Copy, Debug)]
// pub struct Vec3 {
//     pub x: f64,
//     pub y: f64,
//     pub z: f64,
// }
//
//
//
//
// impl Vec3 {
//     pub fn new(x: f64, y: f64, z: f64) -> Self {
//         Self { x, y, z }
//     }
// }
//
