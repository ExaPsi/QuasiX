============
Contributing
============

We welcome contributions to QuasiX! This guide will help you get started.

Ways to Contribute
==================

* **Bug reports**: Open an issue describing the bug
* **Feature requests**: Suggest new features via issues
* **Code contributions**: Submit pull requests
* **Documentation**: Improve docs, examples, tutorials
* **Benchmarks**: Add validation data and comparisons

Development Setup
=================

Prerequisites
-------------

* **Rust**: 1.70 or later
* **Python**: 3.9 or later
* **Git**: For version control

Clone and Setup
---------------

.. code-block:: bash

   # Fork on GitHub, then clone your fork
   git clone git@github.com:YOUR_USERNAME/QuasiX.git
   cd QuasiX

   # Add upstream remote
   git remote add upstream git@github.com:ExaPsi/QuasiX.git

   # Create development environment
   python -m venv .venv
   source .venv/bin/activate
   pip install -e ".[dev]"

Build Rust Extension
--------------------

.. code-block:: bash

   cd quasix
   maturin develop --release
   cd ..

Code Style
==========

Rust
----

* Follow standard Rust conventions
* Run ``cargo fmt`` before committing
* Run ``cargo clippy`` and fix warnings
* Document public functions with ``///`` comments

Example:

.. code-block:: rust

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

Python
------

* Follow PEP 8
* Use type hints
* Run ``ruff check`` and ``ruff format``
* Document with NumPy-style docstrings

Example:

.. code-block:: python

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

Testing
=======

Run Tests
---------

.. code-block:: bash

   # Python tests
   pytest tests/

   # Rust tests
   cargo test --workspace

   # Specific test
   pytest tests/validation/test_g0w0_accuracy.py -v

Add Tests
---------

* Add tests for new features
* Include edge cases
* Test against reference values where possible

Example:

.. code-block:: python

   def test_g0w0_h2o():
       """Test G0W0 for water against reference."""
       mol = gto.M(
           atom='O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587',
           basis='def2-svp'
       )
       mf = scf.RHF(mol).run()

       result = G0W0Driver(mf).kernel()

       assert abs(result.homo_qp - (-12.57)) < 0.05  # Reference value

Pull Request Process
====================

1. **Create a branch**

   .. code-block:: bash

      git checkout -b feature/your-feature

2. **Make changes**

   * Write code
   * Add tests
   * Update documentation

3. **Run checks**

   .. code-block:: bash

      cargo fmt && cargo clippy
      ruff check --fix && ruff format
      pytest tests/

4. **Commit**

   .. code-block:: bash

      git add -A
      git commit -m "feat: Add feature description"

5. **Push and create PR**

   .. code-block:: bash

      git push origin feature/your-feature

   Then open a PR on GitHub.

Commit Message Format
---------------------

Use conventional commits:

* ``feat:`` New feature
* ``fix:`` Bug fix
* ``docs:`` Documentation only
* ``test:`` Adding tests
* ``refactor:`` Code refactoring
* ``perf:`` Performance improvement

Examples:

.. code-block:: text

   feat: Add evGW implementation
   fix: Correct frequency grid normalization
   docs: Update API reference for BSE
   test: Add GW100 validation for first 25 molecules

Code Review
===========

All PRs require review before merging. Reviewers will check:

* Code quality and style
* Test coverage
* Documentation
* Performance impact
* Breaking changes

Reporting Issues
================

Bug Reports
-----------

Include:

1. QuasiX version
2. Python/Rust versions
3. Minimal reproducible example
4. Expected vs actual behavior
5. Full error traceback

Feature Requests
----------------

Include:

1. Use case description
2. Proposed API (if applicable)
3. References to theory/literature

Development Guidelines
======================

Evidence-Based Development
--------------------------

QuasiX follows evidence-based development:

1. **Theory first**: Ground implementation in established quantum chemistry
2. **Test-driven**: Write tests before implementation
3. **Validate against PySCF**: All implementations must match PySCF within tolerance
4. **Document numerical accuracy**: Report deviations and tolerances

Validation Requirements
-----------------------

New implementations must include:

* Comparison with PySCF reference
* Tolerance documentation (typically < 1e-8 Ha)
* GW100 benchmark results (if applicable)

Questions?
==========

* Open an issue with the "question" label
* See existing issues for common questions

License
=======

By contributing, you agree that your contributions will be licensed
under the MIT/Apache-2.0 dual license.
