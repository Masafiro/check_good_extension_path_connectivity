//! theta function representation and isogeny computations
use crate::curves::{Legendre, LegendreProduct, Rosenhain};
use crate::fq::{ExtensionField, Field, Fp2, Fp4};
use deepsize::DeepSizeOf;

// Coordinate compression
// [0, 1, 2, 3, 4, 6, 8, 9, 12, 15] -> [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

const TRANS_INDICES: [[usize; 10]; 15] = [
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    [4, 1, 5, 3, 0, 2, 8, 7, 6, 9],
    [0, 7, 5, 9, 4, 2, 6, 1, 8, 3],
    [4, 7, 2, 9, 0, 5, 8, 1, 6, 3],
    [6, 7, 2, 3, 8, 5, 0, 1, 4, 9],
    [8, 7, 5, 3, 6, 2, 4, 1, 0, 9],
    [6, 1, 5, 9, 8, 2, 0, 7, 4, 3],
    [8, 1, 2, 9, 6, 5, 4, 7, 0, 3],
    [1, 4, 3, 5, 0, 2, 7, 8, 6, 9],
    [2, 3, 6, 7, 5, 8, 0, 1, 4, 9],
    [7, 8, 3, 5, 6, 2, 1, 4, 0, 9],
    [5, 3, 8, 7, 2, 6, 4, 1, 0, 9],
    [7, 4, 9, 2, 0, 5, 1, 8, 6, 3],
    [1, 8, 9, 2, 6, 5, 7, 4, 0, 3],
    [0, 7, 5, 9, 2, 4, 1, 6, 3, 8],
];

// Signs for each transformation (0: +, 1: -, 2: i, 3: -i)
const TRANS_SIGNS: [[usize; 10]; 15] = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 3, 0, 3, 0, 0, 0, 3, 0, 3],
    [0, 0, 0, 2, 0, 0, 0, 0, 0, 2],
    [0, 3, 0, 1, 0, 0, 0, 3, 0, 1],
    [0, 0, 3, 3, 0, 3, 0, 0, 0, 3],
    [0, 3, 3, 2, 0, 3, 0, 3, 0, 2],
    [0, 0, 3, 1, 0, 3, 0, 0, 0, 1],
    [0, 3, 3, 0, 0, 3, 0, 3, 0, 0],
    [0, 3, 0, 3, 0, 0, 0, 3, 0, 1],
    [0, 0, 3, 3, 0, 3, 0, 0, 0, 1],
    [0, 3, 3, 2, 0, 3, 0, 3, 0, 0],
    [0, 3, 3, 2, 0, 3, 0, 3, 0, 0],
    [0, 3, 2, 3, 0, 0, 0, 3, 0, 3],
    [0, 3, 1, 2, 0, 3, 0, 3, 0, 2],
    [0, 0, 0, 2, 0, 0, 0, 0, 0, 2],
];

const ISOGY_INDICES: [[usize; 10]; 4] = [
    [0, 1, 2, 3, 0, 2, 0, 1, 0, 3],
    [1, 0, 3, 2, 1, 3, 1, 0, 1, 2],
    [2, 3, 0, 1, 2, 0, 2, 3, 2, 1],
    [3, 2, 1, 0, 3, 1, 3, 2, 3, 0],
];

const ISOGY_SIGNS: [[i32; 10]; 4] = [
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, -1, -1, 1, 1, -1, -1],
    [1, 1, 1, 1, 1, 1, -1, -1, -1, -1],
    [1, 1, 1, 1, -1, -1, -1, -1, 1, 1],
];

const SPLIT_INDICES: [[usize; 10]; 10] = [
    [9, 2, 1, 8, 6, 5, 4, 7, 3, 0],
    [2, 9, 8, 1, 6, 5, 7, 4, 3, 0],
    [1, 8, 9, 2, 5, 6, 4, 7, 3, 0],
    [8, 1, 2, 9, 5, 6, 7, 4, 3, 0],
    [6, 2, 5, 8, 9, 1, 3, 7, 4, 0],
    [5, 8, 6, 2, 1, 9, 3, 7, 4, 0],
    [4, 7, 1, 8, 3, 5, 9, 2, 6, 0],
    [7, 4, 8, 1, 3, 5, 2, 9, 6, 0],
    [3, 7, 5, 8, 4, 1, 6, 2, 9, 0],
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
];

const SPLIT_SIGNS: [[usize; 10]; 10] = [
    [0, 3, 3, 0, 1, 0, 1, 0, 0, 0],
    [0, 1, 1, 0, 1, 0, 1, 2, 0, 0],
    [0, 1, 1, 0, 1, 2, 1, 0, 0, 0],
    [0, 3, 3, 0, 3, 0, 3, 0, 2, 2],
    [0, 3, 3, 0, 3, 2, 3, 0, 0, 0],
    [0, 1, 1, 0, 3, 0, 3, 0, 0, 0],
    [0, 3, 3, 0, 3, 0, 3, 2, 0, 0],
    [0, 1, 1, 0, 3, 0, 3, 0, 0, 0],
    [0, 1, 1, 2, 1, 0, 1, 0, 0, 2],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
];

#[derive(Copy, Clone, Debug, Eq, PartialEq, DeepSizeOf)]
pub struct Theta<T: ExtensionField> {
    pub null: [T; 10],
}

impl<T: ExtensionField> Theta<T> {
    /// Construct the theta structure from Rosenhain invariants.
    pub fn from_rosenhain(ram: &Rosenhain<T>) -> Self {
        let la = ram.la;
        let mu = ram.mu;
        let nu = ram.nu;

        let mut null = [T::ZERO; 10];

        let inv_prod = (la * nu * (mu - T::ONE) * (la - nu)).inv();

        null[0] = T::ONE;
        null[4] = (mu * (mu - T::ONE) * (la - nu) * inv_prod)
            .sqrt()
            .expect("sqrt failed in from_rosenhain");
        null[6] = ((mu * (nu - T::ONE) * (la - mu)) * la * inv_prod)
            .sqrt()
            .expect("sqrt failed in from_rosenhain");
        null[1] = ((mu * (nu - T::ONE) * (la - T::ONE)) * (la - nu) * inv_prod)
            .sqrt()
            .expect("sqrt failed in from_rosenhain");
        null[2] = ((mu * (la - T::ONE) * (mu - nu)) * nu * inv_prod)
            .sqrt()
            .expect("sqrt failed in from_rosenhain");
        null[5] = null[2] / (nu * null[4]);
        null[8] = null[6] / (la * null[4]);
        let inv_null1 = null[1].inv();
        null[3] = (nu - T::ONE) * (null[4] * null[5]) * inv_null1;
        null[7] = (la - T::ONE) * (null[4] * null[8]) * inv_null1;
        null[9] = (null[0] * null[3] - null[1] * null[2]) / null[8];

        Self { null }
    }

    /// Recover Rosenhain invariants from a squared theta null-point.
    pub fn to_rosenhain(&self) -> Rosenhain<T> {
        let null = &self.null;
        let inv_prod = (null[4] * null[5] * null[8]).inv();
        let la = null[0] * null[5] * null[6] * inv_prod;
        let mu = null[2] * null[4] * null[6] * inv_prod;
        let nu = null[0] * null[2] * null[8] * inv_prod;
        Rosenhain::new(la, mu, nu)
    }

    /// Transform a squared theta null-point by the `m`-th symplectic matrix.
    pub fn transform_nullpoint(&self, m: usize) -> Self {
        let null = &self.null;
        let mut trans_null = [T::ZERO; 10];
        let index_list = &TRANS_INDICES[m];
        let sign_list = &TRANS_SIGNS[m];

        for num in 0..10 {
            let index = index_list[num] as usize;
            let sign = sign_list[num];
            trans_null[index] = match sign {
                0 => null[num],
                1 => T::I * null[num],
                2 => -null[num],
                3 => -T::I * null[num],
                _ => panic!("Invalid sign in transform_nullpoint"),
            };
        }

        Self { null: trans_null }
    }
}

impl<const P: u32, const D: u32> Theta<Fp2<P, D>> {
    /// Find the index of the first zero if the squared theta null-point corresponds to a product of two elliptic curves.
    pub fn find_split_index(&self) -> Option<usize> {
        self.null.iter().position(|&x| x == Fp2::<P, D>::ZERO)
    }

    /// Recover a pair of Legendre invariants from a squared theta null-point.
    pub fn to_legendre_product(&self) -> LegendreProduct<P, D> {
        let clue = self
            .find_split_index()
            .expect("Theta null-point has no zero entries");

        let index_list = &SPLIT_INDICES[clue];
        let sign_list = &SPLIT_SIGNS[clue];

        let mut trans_null = [Fp2::<P, D>::ZERO; 10];
        for num in 0..10 {
            let index = index_list[num] as usize;
            let sign = sign_list[num];
            trans_null[index] = match sign {
                0 => self.null[num],
                1 => Fp2::<P, D>::I * self.null[num],
                2 => -self.null[num],
                3 => -Fp2::<P, D>::I * self.null[num],
                _ => panic!("Invalid sign in to_legendre_product"),
            };
        }

        let la1 = trans_null[0].square() / (trans_null[0].square() - trans_null[1].square());
        let la2 = trans_null[0].square() / (trans_null[0].square() - trans_null[2].square());
        LegendreProduct {
            e1: Legendre::new(la1),
            e2: Legendre::new(la2),
        }
    }

    /// Compute the squared theta null-point of the `m`-th (2,2)-isogenous surface.
    pub fn compute_twoisogeny(&self, m: usize) -> Self {
        let trans_null = self.transform_nullpoint(m).null;
        let sqrt_null = if trans_null[0].is_zero() {
            [
                Fp2::<P, D>::ZERO,
                trans_null[1],
                (trans_null[1] * trans_null[2])
                    .sqrt()
                    .expect("sqrt failure in compute_twoisogeny"),
                (trans_null[1] * trans_null[3])
                    .sqrt()
                    .expect("sqrt failure in compute_twoisogeny"),
            ]
        } else {
            [
                trans_null[0],
                (trans_null[0] * trans_null[1])
                    .sqrt()
                    .expect("sqrt failure in compute_twoisogeny"),
                (trans_null[0] * trans_null[2])
                    .sqrt()
                    .expect("sqrt failure in compute_twoisogeny"),
                (trans_null[0] * trans_null[3])
                    .sqrt()
                    .expect("sqrt failure in compute_twoisogeny"),
            ]
        };

        let mut product_matrix = [[Fp2::<P, D>::ZERO; 4]; 4];
        for k in 0..4 {
            for j in 0..=k {
                let value = sqrt_null[j] * sqrt_null[k];
                product_matrix[j][k] = value;
                product_matrix[k][j] = value;
            }
        }

        let mut image_null = [Fp2::<P, D>::ZERO; 10];
        for num in 0..10 {
            for k in 0..4 {
                let index = ISOGY_INDICES[k][num] as usize;
                image_null[num] += match ISOGY_SIGNS[k][num] {
                    1 => product_matrix[index][k],
                    -1 => -product_matrix[index][k],
                    _ => panic!("Invalid sign in compute_twoisogeny"),
                };
            }
        }

        Self { null: image_null }
    }

    pub fn compute_all_twoisogenies(&self, only_jacobians: bool, is_good: bool) -> Vec<Self> {
        let mut images = Vec::new();
        for m in if is_good { 0..8 } else { 0..15 } {
            let image = self.compute_twoisogeny(m);
            if only_jacobians {
                if let Some(_) = image.find_split_index() {
                    continue;
                }
            }
            images.push(image);
        }
        images
    }

    /// Generate a superspecial genus-2 curves null-point.
    pub fn generate_ssp_theta() -> Self {
        let e = Legendre::<P, D>::generate_ssg_legendre();
        let la = e.la;

        let mut null = [Fp2::<P, D>::ZERO; 10];
        null[0] = Fp2::<P, D>::ONE;
        null[1] = la.sqrt().unwrap();
        null[2] = null[1];
        null[3] = la;
        null[4] = (Fp2::<P, D>::ONE - la).sqrt().unwrap();
        null[5] = null[1] * null[4];
        null[6] = null[4];
        null[7] = null[5];
        null[8] = Fp2::<P, D>::ONE - la;

        let theta = Self { null };

        // Return genus-2 curve, not product of elliptic curves
        theta.compute_all_twoisogenies(true, false)[0]
    }
}

impl<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32> Theta<Fp4<P, D1, D2C0, D2C1>> {
    /// Compute the squared theta null-point of the `m`-th (2,2)-isogenous surface.
    pub fn compute_twoisogeny(&self, m: usize) -> Self {
        let trans_null = self.transform_nullpoint(m).null;
        let sqrt_null = [
            trans_null[0],
            (trans_null[0] * trans_null[1])
                .sqrt()
                .expect("sqrt failure in compute_twoisogeny"),
            (trans_null[0] * trans_null[2])
                .sqrt()
                .expect("sqrt failure in compute_twoisogeny"),
            (trans_null[0] * trans_null[3])
                .sqrt()
                .expect("sqrt failure in compute_twoisogeny"),
        ];

        let mut product_matrix = [[Fp4::<P, D1, D2C0, D2C1>::ZERO; 4]; 4];
        for k in 0..4 {
            for j in 0..=k {
                let value = sqrt_null[j] * sqrt_null[k];
                product_matrix[j][k] = value;
                product_matrix[k][j] = value;
            }
        }

        let mut image_null = [Fp4::<P, D1, D2C0, D2C1>::ZERO; 10];
        for num in 0..10 {
            for k in 0..4 {
                let index = ISOGY_INDICES[k][num] as usize;
                image_null[num] += match ISOGY_SIGNS[k][num] {
                    1 => product_matrix[index][k],
                    -1 => -product_matrix[index][k],
                    _ => panic!("Invalid sign in compute_twoisogeny"),
                };
            }
        }

        Self { null: image_null }
    }

    pub fn compute_all_twoisogenies(&self, is_good: bool) -> Vec<Self> {
        let mut images = Vec::new();
        for m in if is_good { 0..8 } else { 0..15 } {
            let image = self.compute_twoisogeny(m);
            images.push(image);
        }
        images
    }
}
