//! Finite Field Implementations: Fp, Fp2, Fp4.

use deepsize::DeepSizeOf;
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::fmt::{self, Debug, Display};
use std::hash::{Hash, Hasher};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

// =============================================================================
// Scalar Operation Macros
// =============================================================================

/// Scalar Equality Macro (==, !=)
macro_rules! impl_scalar_eq {
    ($Type:ident, $primitive:ty, [$($impl_generics:tt)*], [$($use_generics:tt)*]) => {
        // Type == primitive
        impl<$($impl_generics)*> PartialEq<$primitive> for $Type<$($use_generics)*> {
            fn eq(&self, other: &$primitive) -> bool {
                *self == Self::from(*other)
            }
        }
        // primitive == Type
        impl<$($impl_generics)*> PartialEq<$Type<$($use_generics)*>> for $primitive {
            fn eq(&self, other: &$Type<$($use_generics)*>) -> bool {
                $Type::<$($use_generics)*>::from(*self) == *other
            }
        }
    };
}

/// Scalar Binary Op Macro (+, -, *, /)
macro_rules! impl_scalar_op {
    ($trait:ident, $fn:ident, $Type:ident, $primitive:ty, [$($impl_gen:tt)*], [$($use_gen:tt)*]) => {
        // Type + primitive
        impl<$($impl_gen)*> $trait<$primitive> for $Type<$($use_gen)*> {
            type Output = Self;
            fn $fn(self, rhs: $primitive) -> Self {
                self.$fn(Self::from(rhs))
            }
        }
        // primitive + Type
        impl<$($impl_gen)*> $trait<$Type<$($use_gen)*>> for $primitive {
            type Output = $Type<$($use_gen)*>;
            fn $fn(self, rhs: $Type<$($use_gen)*>) -> $Type<$($use_gen)*> {
                $Type::<$($use_gen)*>::from(self).$fn(rhs)
            }
        }
    };
}

/// Scalar Assign Op Macro (+=, -=, *=, /=)
macro_rules! impl_scalar_assign {
    ($trait:ident, $fn:ident, $Type:ident, $primitive:ty, [$($impl_gen:tt)*], [$($use_gen:tt)*]) => {
        impl<$($impl_gen)*> $trait<$primitive> for $Type<$($use_gen)*> {
            fn $fn(&mut self, rhs: $primitive) {
                self.$fn(Self::from(rhs));
            }
        }
    };
}

/// Master Macro: Apply Everything

macro_rules! impl_scalar_all {
    // ここも [ ] で受け取るように変更
    ($Type:ident, $primitive:ty, [$($impl_gen:tt)*], [$($use_gen:tt)*]) => {
        // 1. Equality
        impl_scalar_eq!($Type, $primitive, [$($impl_gen)*], [$($use_gen)*]);

        // 2. Binary Ops
        impl_scalar_op!(Add, add, $Type, $primitive, [$($impl_gen)*], [$($use_gen)*]);
        impl_scalar_op!(Sub, sub, $Type, $primitive, [$($impl_gen)*], [$($use_gen)*]);
        impl_scalar_op!(Mul, mul, $Type, $primitive, [$($impl_gen)*], [$($use_gen)*]);
        impl_scalar_op!(Div, div, $Type, $primitive, [$($impl_gen)*], [$($use_gen)*]);

        // 3. Assign Ops
        impl_scalar_assign!(AddAssign, add_assign, $Type, $primitive, [$($impl_gen)*], [$($use_gen)*]);
        impl_scalar_assign!(SubAssign, sub_assign, $Type, $primitive, [$($impl_gen)*], [$($use_gen)*]);
        impl_scalar_assign!(MulAssign, mul_assign, $Type, $primitive, [$($impl_gen)*], [$($use_gen)*]);
        impl_scalar_assign!(DivAssign, div_assign, $Type, $primitive, [$($impl_gen)*], [$($use_gen)*]);
    };
}

// =============================================================================
// Trait Definition
// =============================================================================

pub trait Field:
    Sized
    + Copy
    + Clone
    + Debug
    + PartialEq
    + Eq
    + PartialOrd
    + Ord
    + Display
    + Hash
    + From<u32>
    + From<i32>
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
{
    const ZERO: Self;
    const ONE: Self;
    fn is_zero(&self) -> bool;
    fn order() -> u128;
    fn random(rng: &mut impl Rng) -> Self;
    fn inverse(self) -> Option<Self>;
    fn inv(self) -> Self;
    fn sqrt(self) -> Option<Self>;
    fn pow(mut self, mut exp: u64) -> Self {
        let mut result = Self::ONE;
        while exp > 0 {
            if exp & 1 == 1 {
                result *= self;
            }
            self *= self;
            exp >>= 1;
        }
        result
    }
    fn double(&self) -> Self {
        *self + *self
    }
    fn square(&self) -> Self {
        *self * *self
    }
}

pub trait ExtensionField: Field {
    type BaseField: Field;
    const I: Self;
    fn conjugate(self) -> Self;
    fn norm(self) -> Self::BaseField;
}

// =============================================================================
// Fp Definition
// =============================================================================
#[derive(
    Copy, Clone, Debug, Eq, PartialEq, PartialOrd, Ord, Serialize, Deserialize, DeepSizeOf,
)]
pub struct Fp<const P: u32> {
    value: u32,
}

impl<const P: u32> Fp<P> {
    pub const TWO_INV: Fp<P> = (Self { value: 2 }).inv();

    pub const fn new(value: u32) -> Self {
        Self { value: value % P }
    }

    pub const fn value(&self) -> u32 {
        self.value
    }

    pub const fn is_zero(&self) -> bool {
        self.value == 0
    }

    pub const fn add(self, rhs: Self) -> Self {
        let sum = self.value + rhs.value;
        Self {
            value: if sum >= P { sum - P } else { sum },
        }
    }

    pub const fn sub(self, rhs: Self) -> Self {
        Self {
            value: if self.value >= rhs.value {
                self.value - rhs.value
            } else {
                self.value + P - rhs.value
            },
        }
    }

    pub const fn neg(self) -> Self {
        if self.is_zero() {
            self
        } else {
            Self {
                value: P - self.value,
            }
        }
    }

    pub const fn mul(self, rhs: Self) -> Self {
        Self {
            value: ((self.value as u64 * rhs.value as u64) % (P as u64)) as u32,
        }
    }

    pub const fn pow(mut self, mut exp: u64) -> Self {
        let mut result = Self::ONE;
        while exp > 0 {
            if exp & 1 == 1 {
                result = result.mul(self);
            }
            self = self.mul(self);
            exp >>= 1;
        }
        result
    }

    pub const fn inverse(self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        Some(self.pow(P as u64 - 2))
    }

    pub const fn inv(self) -> Self {
        self.inverse().expect("Inverse of zero does not exist")
    }

    pub const fn checked_div(self, rhs: Self) -> Option<Self> {
        match rhs.inverse() {
            Some(inv) => Some(self.mul(inv)),
            None => None,
        }
    }

    pub const fn div(self, rhs: Self) -> Self {
        self.checked_div(rhs).expect("Division by zero")
    }

    /// Legendre symbol (self | P)
    pub const fn legendre_symbol(self) -> i32 {
        let ls = self.pow((P as u64 - 1) / 2).value;
        if ls > 1 { -1 } else { ls as i32 }
    }

    /// Check if the element is a quadratic residue in Fp
    pub const fn is_square(self) -> bool {
        self.legendre_symbol() == 1 || self.value == 0
    }

    /// Find the smallest non-residue modulo p
    pub const fn find_non_residue() -> u32 {
        let mut d = Self { value: 2 };
        while d.is_square() {
            d = d.add(Self::ONE);
        }
        d.value
    }

    /// Square root using Tonelli-Shanks algorithm
    pub const fn sqrt(self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::ZERO);
        }

        // Return None if not a quadratic residue
        if !self.is_square() {
            return None;
        }

        // Simple case: p % 4 == 3
        if P % 4 == 3 {
            let r = self.pow((P as u64 + 1) / 4);
            return Some(r);
        }

        // Find Q and S such that P - 1 = Q * 2^S with Q odd
        let mut q = P - 1;
        let mut s = 0u32;
        while q & 1 == 0 {
            q >>= 1;
            s += 1;
        }

        // Find a quadratic non-residue z
        let z = Self {
            value: Self::find_non_residue(),
        };

        let mut m = s;
        let mut c = z.pow(q as u64);
        let mut t = self.pow(q as u64);
        let mut r = self.pow((q as u64 + 1) / 2);

        while t.value != 1 {
            // Find the least i such that t^(2^i) ≡ 1 (mod P)
            // Note: Since n is a quadratic residue, 0 < i < m is guaranteed
            let mut t2i = t;
            let mut i = 0u32;
            while t2i.value != 1 {
                t2i = t2i.mul(t2i);
                i += 1;
            }

            let b = c.pow(1 << (m - i - 1));
            c = b.mul(b);
            t = t.mul(c);
            r = r.mul(b);
            m = i;
        }

        Some(r)
    }
}

// -----------------------------------------------------------------------------
// Display & Hash
// -----------------------------------------------------------------------------
impl<const P: u32> fmt::Display for Fp<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.value)
    }
}

impl<const P: u32> Hash for Fp<P> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.value.hash(state);
    }
}

// -----------------------------------------------------------------------------
// Type Conversions (From)
// -----------------------------------------------------------------------------
impl<const P: u32> From<u32> for Fp<P> {
    fn from(val: u32) -> Self {
        Self::new(val)
    }
}

impl<const P: u32> From<i32> for Fp<P> {
    fn from(val: i32) -> Self {
        let r = val % (P as i32);
        Self {
            value: if r < 0 {
                (r + P as i32) as u32
            } else {
                r as u32
            },
        }
    }
}

// -----------------------------------------------------------------------------
// Arithmetic Operators
// -----------------------------------------------------------------------------

impl<const P: u32> Add for Fp<P> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        self.add(rhs)
    }
}

impl<const P: u32> AddAssign for Fp<P> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<const P: u32> Sub for Fp<P> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self.sub(rhs)
    }
}

impl<const P: u32> SubAssign for Fp<P> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<const P: u32> Neg for Fp<P> {
    type Output = Self;
    fn neg(self) -> Self {
        self.neg()
    }
}

impl<const P: u32> Mul for Fp<P> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        self.mul(rhs)
    }
}

impl<const P: u32> MulAssign for Fp<P> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<const P: u32> Div for Fp<P> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        self.div(rhs)
    }
}

impl<const P: u32> DivAssign for Fp<P> {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl_scalar_all!(Fp, u32, [const P: u32], [P]);
impl_scalar_all!(Fp, i32, [const P: u32], [P]);

// -----------------------------------------------------------------------------
// Trait Implementation
// -----------------------------------------------------------------------------
impl<const P: u32> Field for Fp<P> {
    const ZERO: Self = Self { value: 0 };
    const ONE: Self = Self { value: 1 };

    fn is_zero(&self) -> bool {
        self.is_zero()
    }

    fn order() -> u128 {
        P as u128
    }

    fn random(rng: &mut impl Rng) -> Self {
        let val: u32 = rng.random_range(0..P);
        Self { value: val }
    }

    fn inverse(self) -> Option<Self> {
        self.inverse()
    }

    fn inv(self) -> Self {
        self.inv()
    }

    fn sqrt(self) -> Option<Self> {
        self.sqrt()
    }
}

// =============================================================================
// Fp2 Definition
// =============================================================================
/// Quadratic extension field over `Fp` defined by the relation `x^2 = D`.
///
/// For primes `p` where `p % 4 == 3` this corresponds to `D = p - 1` (i.e. `x^2 + 1`).
/// Otherwise choose a small quadratic non-residue `D` so that `x^2 - D` is irreducible.
#[derive(
    Copy, Clone, Debug, Eq, PartialEq, PartialOrd, Ord, Serialize, Deserialize, DeepSizeOf,
)]
pub struct Fp2<const P: u32, const D: u32> {
    pub c0: Fp<P>,
    pub c1: Fp<P>,
}

impl<const P: u32, const D: u32> Fp2<P, D> {
    pub const TWO_INV: Self = Self {
        c0: Fp::<P>::TWO_INV,
        c1: Fp::<P>::ZERO,
    };
    pub const D_FP: Fp<P> = Fp::<P>::new(D);

    pub const fn new(c0: Fp<P>, c1: Fp<P>) -> Self {
        Self { c0, c1 }
    }

    pub const fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    pub const fn add(self, rhs: Self) -> Self {
        Self {
            c0: self.c0.add(rhs.c0),
            c1: self.c1.add(rhs.c1),
        }
    }

    pub const fn sub(self, rhs: Self) -> Self {
        Self {
            c0: self.c0.sub(rhs.c0),
            c1: self.c1.sub(rhs.c1),
        }
    }

    pub const fn neg(self) -> Self {
        Self {
            c0: self.c0.neg(),
            c1: self.c1.neg(),
        }
    }

    /// Multiplication in Fp2 defined as Fp[x] / (x^2 - D)
    /// Formula: (a0 + a1*x)(b0 + b1*x) = (a0*b0 + a1*b1*D) + (a0*b1 + a1*b0)*x
    pub const fn mul(self, rhs: Self) -> Self {
        Self {
            c0: self.c0.mul(rhs.c0).add(self.c1.mul(rhs.c1).mul(Self::D_FP)),
            c1: self.c0.mul(rhs.c1).add(self.c1.mul(rhs.c0)),
        }
    }

    pub const fn pow(mut self, mut exp: u64) -> Self {
        let mut result = Self::ONE;
        while exp > 0 {
            if exp & 1 == 1 {
                result = result.mul(self);
            }
            self = self.mul(self);
            exp >>= 1;
        }
        result
    }

    /// Conjugate of the Fp2 element `a - b*x`
    pub const fn conjugate(self) -> Self {
        Self {
            c0: self.c0,
            c1: self.c1.neg(),
        }
    }

    /// Norm of the Fp2 element `a^2 - b^2 * D`
    pub const fn norm(self) -> Fp<P> {
        self.c0
            .mul(self.c0)
            .sub(self.c1.mul(self.c1).mul(Self::D_FP))
    }

    /// Multiplicative inverse using the norm and conjugate
    pub const fn inverse(self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        let norm_inv = match self.norm().inverse() {
            Some(v) => v,
            None => return None,
        };
        Some(Self {
            c0: self.c0.mul(norm_inv),
            c1: self.c1.neg().mul(norm_inv),
        })
    }

    pub const fn inv(self) -> Self {
        self.inverse().expect("Inverse of zero does not exist")
    }

    pub const fn checked_div(self, rhs: Self) -> Option<Self> {
        match rhs.inverse() {
            Some(inv) => Some(self.mul(inv)),
            None => None,
        }
    }

    pub const fn div(self, rhs: Self) -> Self {
        self.checked_div(rhs).expect("Division by zero")
    }

    /// Check if the element is a quadratic residue in Fp2
    pub const fn is_square(self) -> bool {
        self.norm().is_square()
    }

    /// Find a non-residue in Fp2
    pub const fn find_non_residue() -> (u32, u32) {
        let mut i = 0;
        while i < 100 {
            let d = Self {
                c0: Fp::<P> { value: i },
                c1: Fp::<P>::ONE,
            };
            if !d.is_square() {
                return (d.c0.value(), d.c1.value());
            }
            i += 1;
        }
        panic!("Non-residue not found in Fp2 (logic error)");
    }

    /// Square root in Fp2
    /// Find x + y*i such that (x + y*i)^2 = a + b*i
    pub const fn sqrt(self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::ZERO);
        }

        let norm_sqrt = match self.norm().sqrt() {
            Some(v) => v,
            None => return None,
        };

        let a = self.c0;
        let b = self.c1;

        if let Some(x) = a.add(norm_sqrt).mul(Fp::<P>::TWO_INV).sqrt() {
            if !x.is_zero() {
                return Some(Self {
                    c0: x,
                    c1: b.div(x.add(x)),
                });
            }
        }

        if let Some(x) = a.sub(norm_sqrt).mul(Fp::<P>::TWO_INV).sqrt() {
            if !x.is_zero() {
                return Some(Self {
                    c0: x,
                    c1: b.div(x.add(x)),
                });
            }
        }

        Some(Self {
            c0: Fp::ZERO,
            c1: a.div(Self::D_FP).sqrt().unwrap(),
        })
    }
}

// -----------------------------------------------------------------------------
// Display & Hash
// -----------------------------------------------------------------------------
impl<const P: u32, const D: u32> fmt::Display for Fp2<P, D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} + {}*x", self.c0, self.c1)
    }
}

impl<const P: u32, const D: u32> Hash for Fp2<P, D> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.c0.hash(state);
        self.c1.hash(state);
    }
}

// -----------------------------------------------------------------------------
// Type Conversions (From)
// -----------------------------------------------------------------------------
impl<const P: u32, const D: u32> From<u32> for Fp2<P, D> {
    fn from(val: u32) -> Self {
        Self {
            c0: Fp::<P>::from(val),
            c1: Fp::<P>::ZERO,
        }
    }
}

impl<const P: u32, const D: u32> From<i32> for Fp2<P, D> {
    fn from(val: i32) -> Self {
        Self {
            c0: Fp::<P>::from(val),
            c1: Fp::<P>::ZERO,
        }
    }
}

impl<const P: u32, const D: u32> From<(u32, u32)> for Fp2<P, D> {
    fn from(val: (u32, u32)) -> Self {
        Self {
            c0: Fp::<P>::from(val.0),
            c1: Fp::<P>::from(val.1),
        }
    }
}

impl<const P: u32, const D: u32> From<Fp<P>> for Fp2<P, D> {
    fn from(val: Fp<P>) -> Self {
        Self {
            c0: val,
            c1: Fp::<P>::ZERO,
        }
    }
}

// -----------------------------------------------------------------------------
// Arithmetic Operators
// -----------------------------------------------------------------------------
impl<const P: u32, const D: u32> Add for Fp2<P, D> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        self.add(rhs)
    }
}

impl<const P: u32, const D: u32> AddAssign for Fp2<P, D> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<const P: u32, const D: u32> Sub for Fp2<P, D> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self.sub(rhs)
    }
}

impl<const P: u32, const D: u32> SubAssign for Fp2<P, D> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<const P: u32, const D: u32> Neg for Fp2<P, D> {
    type Output = Self;
    fn neg(self) -> Self {
        self.neg()
    }
}

impl<const P: u32, const D: u32> Mul for Fp2<P, D> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        self.mul(rhs)
    }
}

impl<const P: u32, const D: u32> MulAssign for Fp2<P, D> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<const P: u32, const D: u32> Div for Fp2<P, D> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        self.div(rhs)
    }
}

impl<const P: u32, const D: u32> DivAssign for Fp2<P, D> {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl_scalar_all!(Fp2, u32, [const P: u32, const D: u32], [P, D]);
impl_scalar_all!(Fp2, i32, [const P: u32, const D: u32], [P, D]);

// -----------------------------------------------------------------------------
// Trait Implementation
// -----------------------------------------------------------------------------
impl<const P: u32, const D: u32> Field for Fp2<P, D> {
    const ZERO: Self = Self {
        c0: Fp::<P>::ZERO,
        c1: Fp::<P>::ZERO,
    };
    const ONE: Self = Self {
        c0: Fp::<P>::ONE,
        c1: Fp::<P>::ZERO,
    };

    fn is_zero(&self) -> bool {
        self.is_zero()
    }

    fn order() -> u128 {
        let p = P as u128;
        p * p
    }

    fn random(rng: &mut impl Rng) -> Self {
        Self {
            c0: Fp::<P>::random(rng),
            c1: Fp::<P>::random(rng),
        }
    }

    fn inverse(self) -> Option<Self> {
        self.inverse()
    }

    fn inv(self) -> Self {
        self.inv()
    }

    fn sqrt(self) -> Option<Self> {
        self.sqrt()
    }
}

impl<const P: u32, const D: u32> ExtensionField for Fp2<P, D> {
    type BaseField = Fp<P>;

    const I: Self = Fp2::<P, D>::new(Fp::<P>::new(P - 1), Fp::<P>::ZERO)
        .sqrt()
        .unwrap();

    fn conjugate(self) -> Self {
        self.conjugate()
    }

    fn norm(self) -> Self::BaseField {
        self.norm()
    }
}

// =============================================================================
// Fp4 Definition
// =============================================================================
/// Quartic extension field over `Fp2` defined by the relation `y^2 = D2`,
#[derive(
    Copy, Clone, Debug, Eq, PartialEq, PartialOrd, Ord, Serialize, Deserialize, DeepSizeOf,
)]
pub struct Fp4<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> {
    pub c0: Fp2<P, D1>,
    pub c1: Fp2<P, D1>,
}

impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> Fp4<P, D1, D2C0, D2C1> {
    pub const D2_FP2: Fp2<P, D1> = Fp2::<P, D1>::new(Fp::<P>::new(D2C0), Fp::<P>::new(D2C1));

    pub const fn new(c0: Fp2<P, D1>, c1: Fp2<P, D1>) -> Self {
        Self { c0, c1 }
    }

    pub const fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    pub const fn add(self, rhs: Self) -> Self {
        Self {
            c0: self.c0.add(rhs.c0),
            c1: self.c1.add(rhs.c1),
        }
    }

    pub const fn sub(self, rhs: Self) -> Self {
        Self {
            c0: self.c0.sub(rhs.c0),
            c1: self.c1.sub(rhs.c1),
        }
    }

    pub const fn neg(self) -> Self {
        Self {
            c0: self.c0.neg(),
            c1: self.c1.neg(),
        }
    }

    /// Multiplication in Fp4 defined as Fp2[y] / (y^2 - D2)
    /// Formula: (a0 + a1*y)(b0 + b1*y) = (a0*b0 + a1*b1*D2) + (a0*b1 + a1*b0)*y
    pub const fn mul(self, rhs: Self) -> Self {
        Self {
            c0: self
                .c0
                .mul(rhs.c0)
                .add(self.c1.mul(rhs.c1).mul(Self::D2_FP2)),
            c1: self.c0.mul(rhs.c1).add(self.c1.mul(rhs.c0)),
        }
    }

    pub const fn pow(mut self, mut exp: u64) -> Self {
        let mut result = Self::ONE;
        while exp > 0 {
            if exp & 1 == 1 {
                result = result.mul(self);
            }
            self = self.mul(self);
            exp >>= 1;
        }
        result
    }

    /// Conjugate of the Fp4 element `a - b*y`
    pub const fn conjugate(self) -> Self {
        Self {
            c0: self.c0,
            c1: self.c1.neg(),
        }
    }

    /// Norm of the Fp4 element `a^2 - b^2 * D2`
    pub const fn norm(self) -> Fp2<P, D1> {
        self.c0
            .mul(self.c0)
            .sub(self.c1.mul(self.c1).mul(Self::D2_FP2))
    }

    /// Multiplicative inverse using the norm and conjugate
    /// (a + b*y)^-1 = (a - b*y) / (a^2 - b^2*D2)
    pub const fn inverse(self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        let norm_inv = match self.norm().inverse() {
            Some(v) => v,
            None => return None,
        };
        Some(Self {
            c0: self.c0.mul(norm_inv),
            c1: self.c1.neg().mul(norm_inv),
        })
    }

    pub const fn inv(self) -> Self {
        self.inverse().expect("Inverse of zero does not exist")
    }

    pub const fn checked_div(self, rhs: Self) -> Option<Self> {
        match rhs.inverse() {
            Some(inv) => Some(self.mul(inv)),
            None => None,
        }
    }

    pub const fn div(self, rhs: Self) -> Self {
        self.checked_div(rhs).expect("Division by zero")
    }

    /// Check if the element is a quadratic residue in Fp4
    pub const fn sqrt(self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::ZERO);
        }

        let norm_sqrt = match self.norm().sqrt() {
            Some(v) => v,
            None => return None,
        };

        let a = self.c0;
        let b = self.c1;

        if let Some(x) = a.add(norm_sqrt).mul(Fp2::<P, D1>::TWO_INV).sqrt() {
            if !x.is_zero() {
                return Some(Self {
                    c0: x,
                    c1: b.div(x.add(x)),
                });
            }
        }

        if let Some(x) = a.sub(norm_sqrt).mul(Fp2::<P, D1>::TWO_INV).sqrt() {
            if !x.is_zero() {
                return Some(Self {
                    c0: x,
                    c1: b.div(x.add(x)),
                });
            }
        }

        Some(Self {
            c0: Fp2::<P, D1>::ZERO,
            c1: a.div(Self::D2_FP2).sqrt().unwrap(),
        })
    }
}

// -----------------------------------------------------------------------------
// Display & Hash
// -----------------------------------------------------------------------------
impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> fmt::Display
    for Fp4<P, D1, D2C0, D2C1>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} + {}*y", self.c0, self.c1)
    }
}

impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> Hash
    for Fp4<P, D1, D2C0, D2C1>
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.c0.hash(state);
        self.c1.hash(state);
    }
}

// -----------------------------------------------------------------------------
// Type Conversions (From)
// -----------------------------------------------------------------------------
impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> From<u32>
    for Fp4<P, D1, D2C0, D2C1>
{
    fn from(val: u32) -> Self {
        Self {
            c0: Fp2::<P, D1>::from(val),
            c1: Fp2::<P, D1>::ZERO,
        }
    }
}

impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> From<i32>
    for Fp4<P, D1, D2C0, D2C1>
{
    fn from(val: i32) -> Self {
        Self {
            c0: Fp2::<P, D1>::from(val),
            c1: Fp2::<P, D1>::ZERO,
        }
    }
}

impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> From<(u32, u32, u32, u32)>
    for Fp4<P, D1, D2C0, D2C1>
{
    fn from(val: (u32, u32, u32, u32)) -> Self {
        Self {
            c0: Fp2::<P, D1>::new(Fp::<P>::from(val.0), Fp::<P>::from(val.1)),
            c1: Fp2::<P, D1>::new(Fp::<P>::from(val.2), Fp::<P>::from(val.3)),
        }
    }
}

impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> From<Fp<P>>
    for Fp4<P, D1, D2C0, D2C1>
{
    fn from(val: Fp<P>) -> Self {
        Self {
            c0: Fp2::<P, D1>::from(val),
            c1: Fp2::<P, D1>::ZERO,
        }
    }
}

impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> From<Fp2<P, D1>>
    for Fp4<P, D1, D2C0, D2C1>
{
    fn from(val: Fp2<P, D1>) -> Self {
        Self {
            c0: val,
            c1: Fp2::<P, D1>::ZERO,
        }
    }
}

// -----------------------------------------------------------------------------
// Arithmetic Operators
// -----------------------------------------------------------------------------
impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> Add for Fp4<P, D1, D2C0, D2C1> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        self.add(rhs)
    }
}

impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> AddAssign
    for Fp4<P, D1, D2C0, D2C1>
{
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> Sub for Fp4<P, D1, D2C0, D2C1> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self.sub(rhs)
    }
}

impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> SubAssign
    for Fp4<P, D1, D2C0, D2C1>
{
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> Neg for Fp4<P, D1, D2C0, D2C1> {
    type Output = Self;
    fn neg(self) -> Self {
        self.neg()
    }
}

impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> Mul for Fp4<P, D1, D2C0, D2C1> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        self.mul(rhs)
    }
}

impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> MulAssign
    for Fp4<P, D1, D2C0, D2C1>
{
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> Div for Fp4<P, D1, D2C0, D2C1> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        self.div(rhs)
    }
}

impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> DivAssign
    for Fp4<P, D1, D2C0, D2C1>
{
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl_scalar_all!(Fp4, u32, [const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32], [P, D1, D2C0, D2C1]);
impl_scalar_all!(Fp4, i32, [const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32], [P, D1, D2C0, D2C1]);

// -----------------------------------------------------------------------------
// Trait Implementation
// -----------------------------------------------------------------------------
impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> Field
    for Fp4<P, D1, D2C0, D2C1>
{
    const ZERO: Self = Self {
        c0: Fp2::<P, D1>::ZERO,
        c1: Fp2::<P, D1>::ZERO,
    };
    const ONE: Self = Self {
        c0: Fp2::<P, D1>::ONE,
        c1: Fp2::<P, D1>::ZERO,
    };

    fn is_zero(&self) -> bool {
        self.is_zero()
    }

    fn order() -> u128 {
        let p = P as u128;
        p * p * p * p
    }

    fn random(rng: &mut impl Rng) -> Self {
        Self {
            c0: Fp2::<P, D1>::random(rng),
            c1: Fp2::<P, D1>::random(rng),
        }
    }

    fn inverse(self) -> Option<Self> {
        self.inverse()
    }

    fn inv(self) -> Self {
        self.inv()
    }
    fn sqrt(self) -> Option<Self> {
        self.sqrt()
    }
}

impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> ExtensionField
    for Fp4<P, D1, D2C0, D2C1>
{
    type BaseField = Fp2<P, D1>;

    const I: Self = Fp4::<P, D1, D2C0, D2C1>::new(
        Fp2::<P, D1>::new(Fp::<P>::new(P - 1), Fp::<P>::ZERO),
        Fp2::<P, D1>::ZERO,
    )
    .sqrt()
    .unwrap();

    fn conjugate(self) -> Self {
        self.conjugate()
    }

    fn norm(self) -> Self::BaseField {
        self.norm()
    }
}

// =============================================================================
// Tests
// =============================================================================
#[cfg(test)]
mod tests {
    use super::*;

    const P: u32 = 17;
    const D1: u32 = 3;
    const D2C0: u32 = 0;
    const D2C1: u32 = 1;

    #[test]
    fn test_fp_inverse() {
        let a = Fp::<P>::new(5);
        let a_inv = a.inv();
        assert_eq!(a * a_inv, 1);
    }

    #[test]
    fn test_fp_sqrt() {
        let a = Fp::<P>::new(4);
        let a_sq = a * a;
        let a_sqrt = a_sq.sqrt().unwrap();
        assert_eq!(a_sq, a_sqrt * a_sqrt);
    }

    #[test]
    fn test_fp2_inverse() {
        let a = Fp2::<P, D1>::from((5, 1));
        let a_inv = a.inv();
        assert_eq!(a * a_inv, 1);
    }

    #[test]
    fn test_fp2_sqrt() {
        let a = Fp2::<P, D1>::from((5, 1));
        let a_sq = a * a;
        let a_sqrt = a_sq.sqrt().unwrap();
        assert_eq!(a_sq, a_sqrt * a_sqrt);
    }

    #[test]
    fn test_fp4_inverse() {
        let a = Fp4::<P, D1, D2C0, D2C1>::from((2, 1, 3, 4));
        let a_inv = a.inv();
        assert_eq!(a * a_inv, 1);
    }

    #[test]
    fn test_fp4_sqrt() {
        let a = Fp4::<P, D1, D2C0, D2C1>::from((2, 1, 3, 4));
        let a_sq = a * a;
        let a_sqrt = a_sq.sqrt().unwrap();
        assert_eq!(a_sq, a_sqrt * a_sqrt);
    }
}
