# Contributing to QuasiX

We welcome contributions to QuasiX! This guide will help you get started.

## Ways to Contribute

- **Bug reports**: Open an issue describing the bug
- **Feature requests**: Suggest new features via issues
- **Code contributions**: Submit pull requests
- **Documentation**: Improve docs, examples, tutorials
- **Benchmarks**: Add validation data and comparisons

## Development Setup

### Prerequisites

- Rust 1.70+
- Python 3.9+
- Git

### Clone and Setup

```bash
# Fork on GitHub, then clone your fork
git clone git@github.com:YOUR_USERNAME/QuasiX.git
cd QuasiX

# Add upstream remote
git remote add upstream git@github.com:ExaPsi/QuasiX.git

# Create development environment
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

### Build Rust Extension

```bash
cd quasix
maturin develop --release
cd ..
```

## Code Style

### Rust

- Follow standard Rust conventions
- Run `cargo fmt` before committing
- Run `cargo clippy` and fix warnings
- Document public functions with `///` comments

```rust
/// Compute the correlation self-energy using the AC method.
///
/// # Arguments
/// * `energies` - Orbital energies in Hartree
/// * `n_occ` - Number of occupied orbitals
///
/// # Returns
/// Self-energy matrix in Hartree
pub fn compute_sigma_c(energies: &[f64], n_occ: usize) -> Array2<f64> {
    // Implementation
}
```

### Python

- Follow PEP 8
- Use type hints
- Run `ruff check` and `ruff format`
- Document with NumPy-style docstrings

```python
def compute_qp_energies(
    mf: scf.RHF,
    n_freq: int = 32,
) -> np.ndarray:
    """
    Compute quasiparticle energies.

    Parameters
    ----------
    mf : scf.RHF
        PySCF mean-field object.
    n_freq : int, optional
        Number of frequency points, by default 32.

    Returns
    -------
    np.ndarray
        Quasiparticle energies in eV.
    """
    ...
```

## Testing

### Run Tests

```bash
# Python tests
pytest tests/

# Rust tests
cargo test --workspace

# Specific test
pytest tests/validation/test_g0w0_accuracy.py -v
```

### Add Tests

- Add tests for new features
- Include edge cases
- Test against reference values where possible

```python
def test_g0w0_h2o():
    """Test G0W0 for water against reference."""
    mol = gto.M(atom='O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587', basis='def2-svp')
    mf = scf.RHF(mol).run()
    
    result = G0W0Driver(mf).kernel()
    
    assert abs(result.homo_qp - (-12.57)) < 0.05  # Reference value
```

## Pull Request Process

1. **Create a branch**
   ```bash
   git checkout -b feature/your-feature
   ```

2. **Make changes**
   - Write code
   - Add tests
   - Update documentation

3. **Run checks**
   ```bash
   cargo fmt && cargo clippy
   ruff check --fix && ruff format
   pytest tests/
   ```

4. **Commit**
   ```bash
   git add -A
   git commit -m "feat: Add feature description"
   ```

5. **Push and create PR**
   ```bash
   git push origin feature/your-feature
   ```
   Then open a PR on GitHub.

### Commit Message Format

Use conventional commits:

- `feat:` New feature
- `fix:` Bug fix
- `docs:` Documentation only
- `test:` Adding tests
- `refactor:` Code refactoring
- `perf:` Performance improvement

Examples:
```
feat: Add evGW implementation
fix: Correct frequency grid normalization
docs: Update API reference for BSE
test: Add GW100 validation for first 25 molecules
```

## Code Review

All PRs require review before merging. Reviewers will check:

- [ ] Code quality and style
- [ ] Test coverage
- [ ] Documentation
- [ ] Performance impact
- [ ] Breaking changes

## Reporting Issues

### Bug Reports

Include:
1. QuasiX version
2. Python/Rust versions
3. Minimal reproducible example
4. Expected vs actual behavior
5. Full error traceback

### Feature Requests

Include:
1. Use case description
2. Proposed API (if applicable)
3. References to theory/literature

## Questions?

- Open an issue with the "question" label
- Email: quasix@example.com

## License

By contributing, you agree that your contributions will be licensed under the MIT/Apache-2.0 dual license.
