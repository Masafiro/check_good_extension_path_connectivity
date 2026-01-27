//! the number of nodes in isogeny graphs
pub fn number_of_ssg_nodes(p_u32: u32) -> u128 {
    let p = p_u32 as u128;
    let p2 = p * p;
    if matches!(p % 5, 1 | 4) {
        ((p2) * (p2 - 1) * (p2 - 1)) / 2880
    } else {
        1 + ((p2 - 9) * (p2 - 3 * p + 8) * (p2 + 3 * p + 8)) / 2880
    }
}

pub fn number_of_ssp_nodes(p_u32: u32) -> u128 {
    let p = p_u32 as f64;
    let ls_m1: f64 = if p_u32 % 4 == 1 { 1.0 } else { -1.0 };
    let ls_m2: f64 = match p_u32 % 8 {
        1 | 3 => 1.0,
        5 | 7 => -1.0,
        _ => unreachable!(),
    };
    let ls_m3: f64 = if p_u32 % 3 == 1 { 1.0 } else { -1.0 };

    let mut num = ((p - 1.0) * (p + 12.0) * (p + 23.0)) / 2880.0;
    num += ((2.0 * p + 13.0) * (1.0 - ls_m1)) / 96.0;
    num += (1.0 - ls_m2) / 8.0;
    num += ((p + 11.0) * (1.0 - ls_m3)) / 36.0;
    num += ((1.0 - ls_m1) * (1.0 - ls_m3)) / 12.0;
    if p_u32 % 5 == 4 {
        num += 4.0 / 5.0;
    }
    num.round() as u128
}

pub fn number_of_ssp_jacobian_nodes(p_u32: u32) -> u128 {
    let p = p_u32 as f64;
    let ls_m1: f64 = if p_u32 % 4 == 1 { 1.0 } else { -1.0 };
    let ls_m2: f64 = match p_u32 % 8 {
        1 | 3 => 1.0,
        5 | 7 => -1.0,
        _ => unreachable!(),
    };
    let ls_m3: f64 = if p_u32 % 3 == 1 { 1.0 } else { -1.0 };

    let mut num = ((p - 1.0) * (p * p + 25.0 * p + 166.0)) / 2880.0;
    num -= (1.0 - ls_m1) / 32.0;
    num += (1.0 - ls_m2) / 8.0;
    num += (1.0 - ls_m3) / 18.0;
    if p_u32 % 5 == 4 {
        num += 4.0 / 5.0;
    }
    num.round() as u128
}
