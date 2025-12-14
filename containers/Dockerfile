# QuasiX dev/CI image with PySCF and Rust build toolchain
# - Python stack via conda-forge for reliable PySCF wheels
# - Optional Rust+maturin to build QuasiX from source

FROM condaforge/mambaforge:latest

SHELL ["/bin/bash", "-lc"]

# Base scientific Python: pin Python; pin PySCF to a stable minor
RUN mamba install -y -n base -c conda-forge \
    python=3.11 \
    pyscf=2.4.* \
    numpy scipy h5py pandas \
    jupyterlab ipykernel \
    pip && \
    mamba clean -afy

# Tooling to build Rust + Python wheels
RUN apt-get update && \
    apt-get install -y --no-install-recommends curl build-essential pkg-config && \
    rm -rf /var/lib/apt/lists/*

# Rust toolchain + maturin
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y && \
    echo 'source "$HOME/.cargo/env"' >> /etc/profile.d/zz-rust.sh
ENV PATH="/root/.cargo/bin:${PATH}"
RUN pip install --no-cache-dir maturin==1.*

WORKDIR /workspace

# Performance/threading defaults (override per host)
ENV OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    OMP_NUM_THREADS=4 \
    QUASIX_LOG=info \
    PYTHONUNBUFFERED=1

# Optional build from source (uncomment in your fork/CI):
# COPY . /workspace
# RUN maturin build --release -m pyproject.toml && \
#     pip install --no-cache-dir target/wheels/quasix-*.whl

EXPOSE 8888

CMD ["bash", "-lc", "python -c 'import pyscf, sys; print({\"pyscf\": pyscf.__version__}); print(\"QuasiX image ready\")' && tail -f /dev/null"]

