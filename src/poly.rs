//! Polynomial arithmetic over finite fields
use crate::fq::Field;
use rand::Rng;
use std::cmp::max;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Polynomial<F: Field> {
    pub coeffs: Vec<F>,
}

impl<F: Field> Polynomial<F> {
    pub fn new(coeffs: Vec<F>) -> Self {
        let mut poly = Self { coeffs };
        poly.trim();
        poly
    }

    /// Zero polynomial
    pub fn zero() -> Self {
        Self { coeffs: Vec::new() }
    }

    pub fn one() -> Self {
        Self {
            coeffs: vec![F::ONE],
        }
    }

    /// Create a constant polynomial
    pub fn constant(c: F) -> Self {
        if c.is_zero() {
            Self::zero()
        } else {
            Self { coeffs: vec![c] }
        }
    }

    /// Trim leading zero coefficients
    fn trim(&mut self) {
        while let Some(true) = self.coeffs.last().map(|c| c.is_zero()) {
            self.coeffs.pop();
        }
    }

    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    pub fn degree(&self) -> usize {
        if self.is_zero() {
            0 // For convenience, we return 0 here
        } else {
            self.coeffs.len() - 1
        }
    }

    /// Coefficient of x^k
    pub fn coeff(&self, k: usize) -> F {
        if k < self.coeffs.len() {
            self.coeffs[k]
        } else {
            F::ZERO
        }
    }

    pub fn derivative(&self) -> Self {
        if self.coeffs.len() <= 1 {
            return Self::zero();
        }
        let mut res = Vec::with_capacity(self.coeffs.len() - 1);
        for (i, &c) in self.coeffs.iter().enumerate().skip(1) {
            res.push(c * F::from(i as u32));
        }
        Self::new(res)
    }

    /// Returns quotient and remainder (self / divisor, self % divisor)
    pub fn div_rem(&self, divisor: &Self) -> (Self, Self) {
        if divisor.is_zero() {
            panic!("Division by zero polynomial");
        }
        if self.is_zero() {
            return (Self::zero(), Self::zero());
        }

        let div_deg = divisor.degree();
        let self_deg = self.degree();

        if self_deg < div_deg {
            return (Self::zero(), self.clone());
        }

        // Preallocate quotient coefficients; length is degree difference + 1 and init to 0
        let mut quotient_coeffs = vec![F::ZERO; self_deg - div_deg + 1];
        let mut remainder = self.clone();

        let div_lead = divisor.coeffs.last().unwrap();
        let div_lead_inv = div_lead
            .inverse()
            .expect("Leading coeff should be invertible");

        while remainder.degree() >= div_deg && !remainder.is_zero() {
            let rem_deg = remainder.degree();
            let shift = rem_deg - div_deg;
            let scale = *remainder.coeffs.last().unwrap() * div_lead_inv;

            // Record the quotient coefficient: coefficient of x^shift is scale
            quotient_coeffs[shift] = scale;

            // Update remainder in place: remainder -= scale * x^shift * divisor
            for (i, &c) in divisor.coeffs.iter().enumerate() {
                if shift + i < remainder.coeffs.len() {
                    remainder.coeffs[shift + i] -= c * scale;
                }
            }
            remainder.trim();
        }

        let quotient = Self::new(quotient_coeffs);
        (quotient, remainder)
    }

    /// Remainder (self % divisor)
    pub fn rem(&self, divisor: &Self) -> Self {
        self.div_rem(divisor).1
    }

    /// Greatest Common Divisor (GCD) - returns a monic polynomial
    pub fn gcd(a: &Self, b: &Self) -> Self {
        let mut a = a.clone();
        let mut b = b.clone();
        while !b.is_zero() {
            let r = a.rem(&b);
            a = b;
            b = r;
        }
        a.make_monic();
        a
    }

    /// Make the leading coefficient 1
    pub fn make_monic(&mut self) {
        if let Some(&lead) = self.coeffs.last() {
            if lead != F::ONE {
                if let Some(inv) = lead.inverse() {
                    for c in self.coeffs.iter_mut() {
                        *c *= inv;
                    }
                }
            }
        }
    }

    /// Resultant
    pub fn resultant(f: &Self, g: &Self) -> F {
        let n = f.degree();
        let m = g.degree();
        if f.is_zero() || g.is_zero() {
            return F::ZERO;
        }

        let size = n + m;
        let mut mat = vec![vec![F::ZERO; size]; size];

        // Construct the Sylvester matrix
        for i in 0..m {
            for (j, &c) in f.coeffs.iter().enumerate() {
                mat[i][i + j] = c;
            }
        }
        for i in 0..n {
            for (j, &c) in g.coeffs.iter().enumerate() {
                mat[m + i][i + j] = c;
            }
        }

        Self::determinant(mat)
    }

    // Determinant using Gaussian elimination
    fn determinant(mut mat: Vec<Vec<F>>) -> F {
        let n = mat.len();
        let mut det = F::ONE;

        for i in 0..n {
            let mut pivot = i;
            while pivot < n && mat[pivot][i].is_zero() {
                pivot += 1;
            }
            if pivot == n {
                return F::ZERO;
            }

            if pivot != i {
                mat.swap(i, pivot);
                det = -det;
            }

            let pivot_val = mat[i][i];
            det *= pivot_val;
            let inv = pivot_val.inverse().unwrap(); // pivot != 0 so safe

            // normalization
            for j in i..n {
                mat[i][j] *= inv;
            }

            // Elimination
            for k in (i + 1)..n {
                if !mat[k][i].is_zero() {
                    let factor = mat[k][i];
                    for j in i..n {
                        let val = mat[i][j] * factor;
                        mat[k][j] -= val;
                    }
                }
            }
        }
        det
    }

    /// English: Calculate Res(f, f') to check for multiple roots
    pub fn discriminant(&self) -> F {
        let deriv = self.derivative();
        Self::resultant(self, &deriv)
    }

    pub fn pow(&self, mut exp: u128) -> Self {
        let mut base = self.clone();
        let mut result = Polynomial::one();
        while exp > 0 {
            if exp % 2 == 1 {
                result *= &base;
            }
            base = &base * &base;
            exp /= 2;
        }
        result
    }

    /// Calculate (self^exp) % modulus
    pub fn pow_mod(&self, exp: u128, modulus: &Self) -> Self {
        if modulus.is_zero() {
            panic!("Modulo by zero polynomial");
        }
        let mut base = self.rem(modulus);
        let mut result = Self::one();
        let mut e = exp;

        while e > 0 {
            if e & 1 == 1 {
                result = (&result * &base).rem(modulus);
            }
            // base = (base * base) % modulus
            base = (&base * &base).rem(modulus);
            e >>= 1;
        }
        result
    }

    /// Find all distinct roots in F (Cantor-Zassenhaus algorithm)
    pub fn roots(&self, rng: &mut impl Rng) -> Vec<F> {
        match self.degree() {
            0 => return vec![],
            1 => return vec![-self.coeff(0) / self.coeff(1)],
            2 => {
                let (a, b, c) = (self.coeff(2), self.coeff(1), self.coeff(0));
                let d = b * b - F::from(4u32) * a * c;
                if d.is_zero() {
                    return vec![-b / (F::from(2u32) * a)];
                }
                if let Some(d_sqrt) = d.sqrt() {
                    let inv_2a = (F::from(2u32) * a).inv();
                    return vec![(-b + d_sqrt) * inv_2a, (-b - d_sqrt) * inv_2a];
                } else {
                    return vec![];
                }
            }
            _ => {}
        }

        // 1. Take GCD with x^q - x to keep only roots in the base field F
        let x_poly = Self::new(vec![F::ZERO, F::ONE]); // f(x) = x
        let xq_mod_self = x_poly.pow_mod(F::order(), self);
        let xq_minus_x = &xq_mod_self - &x_poly;

        let g = Self::gcd(self, &xq_minus_x);

        // Factor recursively from here
        Self::roots_recursive(&g, rng)
    }

    fn roots_recursive(poly: &Self, rng: &mut impl Rng) -> Vec<F> {
        let deg = poly.degree();
        match deg {
            0 => return vec![],
            1 => return vec![-poly.coeff(0) / poly.coeff(1)],
            2 => {
                let (a, b, c) = (poly.coeff(2), poly.coeff(1), poly.coeff(0));
                let d = b * b - F::from(4u32) * a * c;
                if let Some(d_sqrt) = d.sqrt() {
                    let inv_2a = (F::from(2u32) * a).inv();
                    return vec![(-b + d_sqrt) * inv_2a, (-b - d_sqrt) * inv_2a];
                } else {
                    return vec![];
                }
            }
            _ => {}
        }

        loop {
            // 2. Pick a random polynomial h = x + random constant
            let h = Self::new(vec![F::random(rng), F::ONE]);

            // 3. Compute s(x) = h(x)^((q-1)/2) - 1 mod poly
            //    Roots split by whether this becomes zero or not
            let exp = (F::order() - 1) / 2;
            let h_pow = h.pow_mod(exp, poly);
            let s = &h_pow - &Self::one();

            // 4. Take GCD
            let d = Self::gcd(poly, &s);

            // If it fails, retry with a new random choice
            if d.degree() == 0 || d.degree() == deg {
                continue;
            }

            // poly = d * quotient
            let mut roots = Self::roots_recursive(&d, rng);
            let (quotient, _) = poly.div_rem(&d);
            roots.extend(Self::roots_recursive(&quotient, rng));
            return roots;
        }
    }
}

impl<F: Field> Add for &Polynomial<F> {
    type Output = Polynomial<F>;
    fn add(self, rhs: Self) -> Self::Output {
        let max_len = max(self.coeffs.len(), rhs.coeffs.len());
        let mut res = Vec::with_capacity(max_len);
        for i in 0..max_len {
            let a = self.coeff(i);
            let b = rhs.coeff(i);
            res.push(a + b);
        }
        Polynomial::new(res)
    }
}

impl<F: Field> AddAssign<&Polynomial<F>> for Polynomial<F> {
    fn add_assign(&mut self, rhs: &Polynomial<F>) {
        let max_len = max(self.coeffs.len(), rhs.coeffs.len());
        self.coeffs.resize(max_len, F::ZERO);
        for i in 0..max_len {
            let b = rhs.coeff(i);
            self.coeffs[i] += b;
        }
        self.trim();
    }
}

impl<F: Field> Sub for &Polynomial<F> {
    type Output = Polynomial<F>;
    fn sub(self, rhs: Self) -> Self::Output {
        let max_len = max(self.coeffs.len(), rhs.coeffs.len());
        let mut res = Vec::with_capacity(max_len);
        for i in 0..max_len {
            let a = self.coeff(i);
            let b = rhs.coeff(i);
            res.push(a - b);
        }
        Polynomial::new(res)
    }
}

impl<F: Field> SubAssign<&Polynomial<F>> for Polynomial<F> {
    fn sub_assign(&mut self, rhs: &Polynomial<F>) {
        let max_len = max(self.coeffs.len(), rhs.coeffs.len());
        self.coeffs.resize(max_len, F::ZERO);
        for i in 0..max_len {
            let b = rhs.coeff(i);
            self.coeffs[i] -= b;
        }
        self.trim();
    }
}

impl<F: Field> Neg for &Polynomial<F> {
    type Output = Polynomial<F>;
    fn neg(self) -> Self::Output {
        let res: Vec<F> = self.coeffs.iter().map(|&c| -c).collect();
        Polynomial::new(res)
    }
}

impl<F: Field> Mul for &Polynomial<F> {
    type Output = Polynomial<F>;
    fn mul(self, rhs: Self) -> Self::Output {
        if self.is_zero() || rhs.is_zero() {
            return Polynomial::zero();
        }
        let len = self.coeffs.len() + rhs.coeffs.len() - 1;
        let mut res = vec![F::ZERO; len];
        for (i, &c1) in self.coeffs.iter().enumerate() {
            for (j, &c2) in rhs.coeffs.iter().enumerate() {
                res[i + j] += c1 * c2;
            }
        }
        Polynomial::new(res)
    }
}

impl<F: Field> MulAssign<&Polynomial<F>> for Polynomial<F> {
    fn mul_assign(&mut self, rhs: &Polynomial<F>) {
        *self = &*self * rhs;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fq::Fp; // Import existing Fp
    use rand::SeedableRng;

    // Define a polynomial type over P=7 for tests
    type F = Fp<7>;
    type Poly = Polynomial<F>;

    // Helper to construct Fp from integer
    fn n(v: u32) -> F {
        F::new(v)
    }

    #[test]
    fn test_normalize() {
        // [1, 2, 0, 0] -> [1, 2] (trim trailing zeros)
        let p = Poly::new(vec![n(1), n(2), n(0), n(0)]);
        assert_eq!(p.coeffs, vec![n(1), n(2)]);
        assert_eq!(p.degree(), 1);

        // [0, 0] -> [] (zero polynomial)
        let zero = Poly::new(vec![n(0), n(0)]);
        assert!(zero.is_zero());
        assert!(zero.coeffs.is_empty());
    }

    #[test]
    fn test_arithmetic() {
        // p1 = x + 1 ([1, 1])
        let p1 = Poly::new(vec![n(1), n(1)]);
        // p2 = x + 6 (equivalently x - 1) ([6, 1])
        let p2 = Poly::new(vec![n(6), n(1)]);

        // Addition: (x + 1) + (x + 6) = 2x + 7 = 2x ([0, 2]) in F7
        let sum = &p1 + &p2;
        assert_eq!(sum.coeffs, vec![n(0), n(2)]);

        // Multiplication: (x + 1)(x - 1) = x^2 - 1 = x^2 + 6 ([6, 0, 1]) in F7
        let prod = &p1 * &p2;
        assert_eq!(prod.coeffs, vec![n(6), n(0), n(1)]);
    }

    #[test]
    fn test_derivative() {
        // f = 3x^2 + 2x + 5 ([5, 2, 3])
        let f = Poly::new(vec![n(5), n(2), n(3)]);

        // f' = 6x + 2 ([2, 6])
        let diff = f.derivative();
        assert_eq!(diff.coeffs, vec![n(2), n(6)]);
    }

    #[test]
    fn test_gcd() {
        // A = (x-1)(x-2) = x^2 - 3x + 2 = x^2 + 4x + 2
        let a = Poly::new(vec![n(2), n(4), n(1)]);

        // B = (x-1)(x-3) = x^2 - 4x + 3 = x^2 + 3x + 3
        let b = Poly::new(vec![n(3), n(3), n(1)]);

        // GCD should be x-1 = x+6
        let gcd = Poly::gcd(&a, &b);
        assert_eq!(gcd.coeffs, vec![n(6), n(1)]);
    }

    #[test]
    fn test_roots() {
        // Fixed seed RNG for reproducibility
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        // f(x) = (x - 1)(x - 2)(x - 3) = x^3 - 6x^2 + 11x - 6 = x^3 + x^2 + 4x + 1 mod 7
        let p = Poly::new(vec![n(1), n(4), n(1), n(1)]);

        let mut roots = p.roots(&mut rng);
        roots.sort_by_key(|a| a.value()); // Sort results for comparison

        // Expected roots: 1, 2, 3
        assert_eq!(roots, vec![n(1), n(2), n(3)]);
    }

    #[test]
    fn test_roots_with_no_solution() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        // g(x) = x^2 + 1 (no solution in P=7; since 3^2=2 and 4^2=2, -1(=6) never occurs)
        let p = Poly::new(vec![n(1), n(0), n(1)]);

        let roots = p.roots(&mut rng);
        assert!(roots.is_empty());
    }
}
