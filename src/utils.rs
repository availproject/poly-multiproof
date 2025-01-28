//! Crate-wide utility functions.

/// Finds the smallest power of 2 greater than or equal to `a`.
pub fn smallest_power_of_2_greater_than(a: usize) -> usize {
    if a == 0 {
        return 1;
    }

    let mut power = 1;
    while power <= a {
        power <<= 1;
    }

    power
}
