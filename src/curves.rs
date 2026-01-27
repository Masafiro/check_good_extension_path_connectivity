//! genus-2 curves and elliptic curves over Fp2 and Fp4
use crate::{
    fq::{Field, Fp, Fp2, Fp4},
    poly::Polynomial,
};
use core::panic;
use deepsize::DeepSizeOf;
use serde::{Deserialize, Serialize};
use std::cmp::{max, min};

const G2: [[usize; 6]; 15] = [
    [0, 1, 2, 3, 4, 5],
    [2, 0, 1, 3, 4, 5],
    [3, 0, 1, 2, 4, 5],
    [4, 0, 1, 2, 3, 5],
    [5, 0, 1, 2, 3, 4],
    [4, 1, 0, 2, 3, 5],
    [5, 1, 0, 2, 3, 4],
    [4, 2, 0, 1, 3, 5],
    [5, 2, 0, 1, 3, 4],
    [3, 0, 4, 1, 2, 5],
    [4, 0, 3, 1, 2, 5],
    [5, 0, 3, 1, 2, 4],
    [5, 1, 3, 0, 2, 4],
    [2, 3, 4, 0, 1, 5],
    [5, 0, 4, 1, 2, 3],
];

const G4: [[usize; 6]; 10] = [
    [0, 1, 2, 3, 4, 5],
    [3, 0, 1, 2, 4, 5],
    [4, 0, 1, 2, 3, 5],
    [5, 0, 1, 2, 3, 4],
    [2, 3, 0, 1, 4, 5],
    [4, 2, 0, 1, 3, 5],
    [5, 2, 0, 1, 3, 4],
    [3, 4, 0, 1, 2, 5],
    [5, 3, 0, 1, 2, 4],
    [4, 5, 0, 1, 2, 3],
];

const G6: [[usize; 6]; 60] = [
    [0, 1, 2, 3, 4, 5],
    [1, 0, 2, 3, 4, 5],
    [2, 0, 1, 3, 4, 5],
    [3, 0, 1, 2, 4, 5],
    [4, 0, 1, 2, 3, 5],
    [5, 0, 1, 2, 3, 4],
    [0, 2, 1, 3, 4, 5],
    [1, 2, 0, 3, 4, 5],
    [2, 1, 0, 3, 4, 5],
    [3, 1, 0, 2, 4, 5],
    [4, 1, 0, 2, 3, 5],
    [5, 1, 0, 2, 3, 4],
    [0, 3, 1, 2, 4, 5],
    [1, 3, 0, 2, 4, 5],
    [2, 3, 0, 1, 4, 5],
    [3, 2, 0, 1, 4, 5],
    [4, 2, 0, 1, 3, 5],
    [5, 2, 0, 1, 3, 4],
    [0, 4, 1, 2, 3, 5],
    [1, 4, 0, 2, 3, 5],
    [2, 4, 0, 1, 3, 5],
    [3, 4, 0, 1, 2, 5],
    [4, 3, 0, 1, 2, 5],
    [5, 3, 0, 1, 2, 4],
    [0, 5, 1, 2, 3, 4],
    [1, 5, 0, 2, 3, 4],
    [2, 5, 0, 1, 3, 4],
    [3, 5, 0, 1, 2, 4],
    [4, 5, 0, 1, 2, 3],
    [5, 4, 0, 1, 2, 3],
    [0, 1, 3, 2, 4, 5],
    [1, 0, 3, 2, 4, 5],
    [2, 0, 3, 1, 4, 5],
    [3, 0, 2, 1, 4, 5],
    [4, 0, 2, 1, 3, 5],
    [5, 0, 2, 1, 3, 4],
    [0, 2, 3, 1, 4, 5],
    [1, 2, 3, 0, 4, 5],
    [2, 1, 3, 0, 4, 5],
    [0, 3, 2, 1, 4, 5],
    [1, 3, 2, 0, 4, 5],
    [2, 3, 1, 0, 4, 5],
    [0, 4, 2, 1, 3, 5],
    [1, 4, 2, 0, 3, 5],
    [2, 4, 1, 0, 3, 5],
    [0, 5, 2, 1, 3, 4],
    [1, 5, 2, 0, 3, 4],
    [2, 5, 1, 0, 3, 4],
    [0, 1, 4, 2, 3, 5],
    [1, 0, 4, 2, 3, 5],
    [2, 0, 4, 1, 3, 5],
    [3, 0, 4, 1, 2, 5],
    [0, 2, 4, 1, 3, 5],
    [1, 2, 4, 0, 3, 5],
    [2, 1, 4, 0, 3, 5],
    [3, 1, 4, 0, 2, 5],
    [0, 3, 4, 1, 2, 5],
    [1, 3, 4, 0, 2, 5],
    [2, 3, 4, 0, 1, 5],
    [3, 2, 4, 0, 1, 5],
];

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, DeepSizeOf)]
pub struct IgusaInvariants<T: Field> {
    pub i1: T,
    pub i2: T,
    pub i3: T,
}

/// A Rosenhain form of a genus 2 curve: y^2 = x (x - 1) (x - λ) (x - μ) (x - ν)
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize, DeepSizeOf)]
pub struct Rosenhain<T: Field> {
    pub la: T,
    pub mu: T,
    pub nu: T,
}

impl<T: Field> Rosenhain<T> {
    pub fn new(la: T, mu: T, nu: T) -> Self {
        Self { la, mu, nu }
    }

    pub fn igusa_invariants(&self) -> IgusaInvariants<T> {
        // Find a field element 'a' different from la, mu, nu
        let mut a = T::ONE.double();
        while a == self.la || a == self.mu || a == self.nu {
            a = a + T::ONE;
        }

        let rs = [
            T::ZERO,
            T::ONE,
            (T::ONE - a).inv(),
            self.la * (self.la - a).inv(),
            self.mu * (self.mu - a).inv(),
            self.nu * (self.nu - a).inv(),
        ];

        let mut i2 = T::ZERO;
        for g in G2.iter() {
            let t = (rs[g[0]] - rs[g[1]]) * (rs[g[2]] - rs[g[3]]) * (rs[g[4]] - rs[g[5]]);
            i2 += t * t;
        }

        let mut i4 = T::ZERO;
        for g in G4.iter() {
            let t = (rs[g[0]] - rs[g[1]])
                * (rs[g[1]] - rs[g[2]])
                * (rs[g[2]] - rs[g[0]])
                * (rs[g[3]] - rs[g[4]])
                * (rs[g[4]] - rs[g[5]])
                * (rs[g[5]] - rs[g[3]]);
            i4 += t * t;
        }

        let mut i6 = T::ZERO;
        for g in G6.iter() {
            let t = (rs[g[0]] - rs[g[1]])
                * (rs[g[1]] - rs[g[2]])
                * (rs[g[2]] - rs[g[0]])
                * (rs[g[3]] - rs[g[4]])
                * (rs[g[4]] - rs[g[5]])
                * (rs[g[5]] - rs[g[3]])
                * (rs[g[0]] - rs[g[3]])
                * (rs[g[1]] - rs[g[4]])
                * (rs[g[2]] - rs[g[5]]);
            i6 += t * t;
        }

        let mut i10 = T::ONE;
        for i in 0..6 {
            for j in (i + 1)..6 {
                let diff = rs[i] - rs[j];
                i10 *= diff * diff;
            }
        }

        let j2 = i2 / T::from(8u32);
        let j4 = (T::from(4u32) * j2 * j2 - i4) / T::from(96u32);
        let j6 = (T::from(8u32) * j2.pow(3) - T::from(160u32) * j2 * j4 - i6) / T::from(576u32);
        let j10 = i10 / T::from(4096u32);

        if !j2.is_zero() {
            IgusaInvariants {
                i1: j2.pow(5) / j10,
                i2: j2.pow(3) * j4 / j10,
                i3: j2.pow(2) * j6 / j10,
            }
        } else if !j4.is_zero() {
            IgusaInvariants {
                i1: T::ZERO,
                i2: j4.pow(5) / j10.pow(2),
                i3: j4 * j6 / j10,
            }
        } else {
            IgusaInvariants {
                i1: T::ZERO,
                i2: T::ZERO,
                i3: j6.pow(5) / j10.pow(3),
            }
        }
    }
}

impl<const P: u32, const D: u32> Legendre<P, D> {
    pub fn is_supersingular(&self) -> bool {
        let m = (P - 1) / 2;
        let mut binom = Fp::<P>::ONE;
        let mut la_pow = Fp2::<P, D>::ONE;
        let mut sum = Fp2::<P, D>::ZERO;
        for k in 0..=m {
            let binom_square = Fp2::<P, D>::from(binom * binom);
            sum += binom_square * la_pow;
            binom *= Fp::from(m - k) / Fp::from(k + 1);
            la_pow *= self.la;
        }
        sum.is_zero()
    }

    /// Generate a supersingular elliptic curve over Fp^2
    pub fn generate_ssg_legendre() -> Legendre<P, D> {
        assert!(P > 5, "Prime P must be greater than 5");

        if P % 4 == 3 {
            let la = Fp2::<P, D>::new(Fp::from(P - 1), Fp::ZERO);
            return Legendre::<P, D>::new(la);
        }

        for c0 in 0..P {
            for c1 in 0..P {
                let la = Fp2::<P, D>::from((c0, c1));
                if la.is_zero() || la == Fp2::<P, D>::ONE {
                    continue;
                }
                let e = Legendre::<P, D>::new(la);
                if e.is_supersingular() {
                    return e;
                }
            }
        }
        panic!("Failed to find a supersingular elliptic curve over Fp2");
    }
}

impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32>
    Rosenhain<Fp4<P, D1, D2C0, D2C1>>
{
    /// Generate a supersingular genus-2 curve over Fp^4
    pub fn generate_ssg_rosenhain() -> Self {
        assert!(P > 5, "Prime P must be greater than 5");

        let mut rng = rand::rng();

        let f = if P % 5 == 2 || P % 5 == 3 {
            // x^5 - 1
            Polynomial::<Fp4<P, D1, D2C0, D2C1>>::new(vec![
                -Fp4::<P, D1, D2C0, D2C1>::ONE,
                Fp4::<P, D1, D2C0, D2C1>::ZERO,
                Fp4::<P, D1, D2C0, D2C1>::ZERO,
                Fp4::<P, D1, D2C0, D2C1>::ZERO,
                Fp4::<P, D1, D2C0, D2C1>::ZERO,
                Fp4::<P, D1, D2C0, D2C1>::ONE,
            ])
        } else {
            let exp = ((P - 1) / 2) as u128;

            let f_fp = loop {
                let poly = Polynomial::<Fp<P>>::new(vec![
                    Fp::<P>::random(&mut rng),
                    Fp::<P>::random(&mut rng),
                    Fp::<P>::random(&mut rng),
                    Fp::<P>::random(&mut rng),
                    Fp::<P>::ZERO,
                    Fp::<P>::ONE,
                ]);

                if poly.discriminant().is_zero() {
                    continue;
                }

                let g = poly.pow(exp);

                let a = g.coeff((P - 1) as usize);
                let b = g.coeff((2 * P - 1) as usize);
                let c = g.coeff((P - 2) as usize);
                let d = g.coeff((2 * P - 2) as usize);

                // Non-superspecial condition (Cartier-Manin matrix is not zero)
                if a.is_zero() && b.is_zero() && c.is_zero() && d.is_zero() {
                    continue;
                }
                // Supersingular condition (Cartier-Manin matrix is nilpotent)
                if !(a * d - b * c).is_zero() || !(a + d).is_zero() {
                    continue;
                }
                break poly;
            };
            let coeffs_fp_4 = f_fp
                .coeffs
                .iter()
                .map(|&c| Fp4::<P, D1, D2C0, D2C1>::from(c))
                .collect();
            Polynomial::<Fp4<P, D1, D2C0, D2C1>>::new(coeffs_fp_4)
        };
        let roots = f.roots(&mut rng);
        assert!(roots.len() == 5, "Failed to find 5 distinct roots");

        let debom_inv = (roots[1] - roots[0]).inv();
        let la = (roots[2] - roots[0]) * debom_inv;
        let mu = (roots[3] - roots[0]) * debom_inv;
        let nu = (roots[4] - roots[0]) * debom_inv;

        Self { la, mu, nu }
    }
}

/// A Legendre form of an elliptic curve: y^2 = x (x - 1) (x - λ)
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct Legendre<const P: u32, const D: u32> {
    pub la: Fp2<P, D>,
}

impl<const P: u32, const D: u32> Legendre<P, D> {
    pub fn new(la: Fp2<P, D>) -> Self {
        Self { la }
    }

    pub fn j_invariant(&self) -> Fp2<P, D> {
        let la = self.la;
        let one = Fp2::<P, D>::ONE;
        (256 * (la * la - la + one).pow(3)) / (la * la * (la - one) * (la - one))
    }
}

/// Two legendre forms of elliptic curves.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct LegendreProduct<const P: u32, const D: u32> {
    pub e1: Legendre<P, D>,
    pub e2: Legendre<P, D>,
}

impl<const P: u32, const D: u32> LegendreProduct<P, D> {
    pub fn new(e1: Legendre<P, D>, e2: Legendre<P, D>) -> Self {
        Self { e1, e2 }
    }

    pub fn j_invariants(&self) -> (Fp2<P, D>, Fp2<P, D>) {
        let j1 = self.e1.j_invariant();
        let j2 = self.e2.j_invariant();
        (min(j1, j2), max(j1, j2))
    }
}
