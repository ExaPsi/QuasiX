use ndarray::Array2;

fn create_spd_matrix(n: usize, condition: f64) -> Array2<f64> {
    let mut q = Array2::zeros((n, n));
    for i in 0..n {
        for j in 0..n {
            q[[i, j]] = if i == j { 1.0 } else { 0.0 };
        }
    }

    let mut d = Array2::zeros((n, n));
    for i in 0..n {
        d[[i, i]] = 1.0 + (condition - 1.0) * (i as f64) / ((n - 1) as f64);
    }

    let temp = q.dot(&d);
    let mut matrix = temp.dot(&q.t());

    for i in 0..n {
        for j in i + 1..n {
            let avg = (matrix[[i, j]] + matrix[[j, i]]) / 2.0;
            matrix[[i, j]] = avg;
            matrix[[j, i]] = avg;
        }
    }

    // Print diagonal elements
    for i in 0..n {
        println!("diag[{}] = {:.2e}", i, matrix[[i, i]]);
    }

    matrix
}

fn main() {
    let matrix = create_spd_matrix(5, 1e8);

    // Compute condition from diagonal
    let mut min_diag = f64::INFINITY;
    let mut max_diag = f64::NEG_INFINITY;
    for i in 0..5 {
        min_diag = min_diag.min(matrix[[i, i]]);
        max_diag = max_diag.max(matrix[[i, i]]);
    }

    println!("\nCondition from diagonal: {:.2e}", max_diag / min_diag);
}
