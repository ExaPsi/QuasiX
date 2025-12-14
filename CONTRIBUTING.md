# Contributing to QuasiX

Thank you for your interest in contributing to QuasiX! This document provides guidelines for contributing to the project.

## 🎯 Development Philosophy

1. **Performance First**: Optimize hot paths, use SIMD where beneficial
2. **Correctness**: Validate against PySCF and published benchmarks
3. **Clean Code**: Follow Rust and Python idioms
4. **Documentation**: Document theory, implementation, and usage
5. **Testing**: Comprehensive tests with clear acceptance criteria

## 🏗️ Project Structure

```
QuasiX/
├── quasix_core/       # Rust computational kernel
│   └── src/
│       ├── df/        # Density fitting
│       ├── gw/        # GW implementation
│       ├── bse/       # BSE solver
│       └── ...
├── quasix/            # Python bindings
│   ├── src/           # Rust PyO3 code
│   └── quasix/        # Python modules
├── tests/             # Test infrastructure
│   ├── verification/  # Story verification
│   └── logs/         # Test outputs (gitignored)
├── docs/             # Documentation
└── examples/         # Usage examples
```

## 🔧 Development Setup

### Prerequisites
- Rust 1.75+ (install via [rustup](https://rustup.rs/))
- Python 3.10+
- Git

### Setup Steps

```bash
# Clone the repository
git clone https://github.com/quasix/quasix.git
cd QuasiX

# Set up Python virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install maturin numpy scipy pyscf pytest black ruff mypy

# Build the project
cd quasix
maturin develop --release
cd ..
```

## 📝 Code Style

### Rust Code
- Follow standard Rust conventions
- Use `cargo fmt` before committing
- Run `cargo clippy` and address warnings
- Document public APIs with doc comments

```rust
/// Compute the exchange self-energy
///
/// # Arguments
/// * `iaP` - Transition density fitting tensor
/// 
/// # Returns
/// Exchange self-energy diagonal elements
pub fn compute_exchange_selfenergy(&self, iaP: &Array2<f64>) -> Result<Array1<f64>> {
    // Implementation
}
```

### Python Code
- Follow PEP 8 conventions
- Use type hints for function signatures
- Format with `black`
- Lint with `ruff`
- Type check with `mypy`

```python
def compute_gw_correction(
    mo_energy: np.ndarray,
    self_energy: np.ndarray,
    z_factor: float
) -> np.ndarray:
    """Compute GW quasiparticle correction.
    
    Parameters
    ----------
    mo_energy : array_like
        Molecular orbital energies
    self_energy : array_like
        Self-energy matrix elements
    z_factor : float
        Quasiparticle weight
        
    Returns
    -------
    array_like
        Quasiparticle energies
    """
    # Implementation
```

## 🧪 Testing Guidelines

### Test Categories

1. **Unit Tests**: Test individual functions
2. **Integration Tests**: Test module interactions
3. **Verification Tests**: Validate acceptance criteria
4. **Benchmark Tests**: Compare with reference data

### Writing Tests

#### Rust Tests
```rust
#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_exchange_selfenergy() {
        let iaP = Array2::zeros((10, 20));
        let result = compute_exchange_selfenergy(&iaP).unwrap();
        assert_eq!(result.len(), 10);
    }
}
```

#### Python Tests
```python
import pytest
import numpy as np
from quasix import compute_gw_correction

def test_gw_correction():
    mo_energy = np.array([1.0, 2.0, 3.0])
    self_energy = np.array([0.1, 0.2, 0.3])
    z_factor = 0.8
    
    result = compute_gw_correction(mo_energy, self_energy, z_factor)
    
    assert result.shape == mo_energy.shape
    np.testing.assert_allclose(result, expected, rtol=1e-6)
```

### Verification Scripts

Create verification scripts in `tests/verification/` following the naming pattern `verify_s{sprint}-{story}.sh`:

```bash
#!/bin/bash
# verify_s2-1.sh - Verification for S2-1: DF tensor construction

set -e

# Set up paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
LOG_DIR="$ROOT_DIR/tests/logs"

mkdir -p "$LOG_DIR"

echo "=== S2-1 Verification ==="
# Test implementation...
```

## 🌟 Contribution Workflow

1. **Check Issues**: Look for open issues or create one for discussion
2. **Fork & Branch**: Create a feature branch from `main`
3. **Implement**: Follow the coding standards
4. **Test**: Add tests and run existing test suite
5. **Document**: Update relevant documentation
6. **PR**: Submit pull request with clear description

### Commit Messages

Follow conventional commits:
```
feat: Add contour deformation integration
fix: Correct sign in exchange self-energy
docs: Update GW theory documentation
test: Add BSE kernel validation tests
perf: Optimize DF tensor contraction with BLAS
```

## 📊 Performance Considerations

When optimizing performance:

1. **Profile First**: Use `cargo flamegraph` or `py-spy`
2. **Benchmark**: Create benchmarks before optimizing
3. **SIMD**: Use explicit SIMD for hot loops when beneficial
4. **Parallelization**: Use Rayon for Rust, consider GIL release for Python
5. **Memory**: Minimize allocations in hot paths

## 📚 Documentation

### Where to Document

- **API Docs**: In code with doc comments
- **Theory**: `docs/theory.md`
- **Implementation**: `docs/stories/`
- **Examples**: `examples/` with working code
- **User Guide**: In README and dedicated guides

### Documentation Checklist

- [ ] Public APIs have doc comments
- [ ] Complex algorithms have explanation
- [ ] Examples demonstrate usage
- [ ] Theory references are cited

## 🐛 Reporting Issues

When reporting issues, please include:

1. **Description**: Clear problem statement
2. **Reproduction**: Minimal code to reproduce
3. **Environment**: OS, Python version, Rust version
4. **Expected**: What should happen
5. **Actual**: What actually happens
6. **Logs**: Relevant error messages

## 🤝 Code Review Process

Pull requests will be reviewed for:

1. **Correctness**: Tests pass, logic is sound
2. **Performance**: No unnecessary overhead
3. **Style**: Follows project conventions
4. **Documentation**: APIs are documented
5. **Tests**: Adequate test coverage

## 📮 Communication

- **Issues**: GitHub issues for bugs and features
- **Discussions**: GitHub discussions for questions
- **Email**: dev@quasix.org for private matters

## 📜 License

By contributing to QuasiX, you agree that your contributions will be licensed under the same terms as the project (MIT/Apache 2.0 dual license).

## 🙏 Recognition

Contributors will be recognized in:
- AUTHORS.md file
- Release notes
- Project documentation

Thank you for contributing to QuasiX!