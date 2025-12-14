"""
QuasiX I/O module - HDF5 support for DF tensors

This module provides high-level Python interfaces to the Rust-based HDF5 I/O
for density-fitting tensors following the S2.4 production schema.
"""

import numpy as np
from typing import Optional

# Try to import Rust HDF5 backend
HAS_RUST_HDF5 = False
DFTensorsHDF5 = None
DFMetadata = None

try:
    # Import from quasix_core (the compiled Rust extension)
    from quasix.quasix import DFTensorsHDF5, DFMetadata
    HAS_RUST_HDF5 = True
except ImportError as e:
    import warnings
    warnings.warn(
        f"Rust HDF5 backend not available: {e}\n"
        "DF tensor HDF5 I/O will fall back to pure Python implementation."
    )

__all__ = ['save_df_tensors_hdf5', 'load_df_tensors_hdf5', 'HAS_RUST_HDF5',
           'DFTensorsHDF5', 'DFMetadata']


def save_df_tensors_hdf5(df_tensor, filename: str) -> None:
    """
    Save DFTensor to HDF5 using Rust implementation

    Args:
        df_tensor: DFTensor object from quasix.df_tensor
        filename: Output HDF5 file path

    Raises:
        ImportError: If Rust HDF5 backend is not available
        IOError: If HDF5 write fails
    """
    if not HAS_RUST_HDF5:
        raise ImportError(
            "Rust HDF5 backend not available. "
            "Please rebuild with: maturin develop --release --features hdf5_support"
        )

    # Validate required attributes
    required_attrs = ['naux', 'nao', 'nocc', 'nvir', 'auxbasis']
    for attr in required_attrs:
        if not hasattr(df_tensor, attr) or getattr(df_tensor, attr) is None:
            raise ValueError(f"DFTensor missing required attribute: {attr}")

    # Get tensors
    metric_2c = df_tensor.get_metric()
    cderi_3c = df_tensor.get_cholesky_vectors()
    ia_p = df_tensor.get_transition_tensor()
    ij_p = df_tensor.get_occupied_tensor()

    # Validate tensors
    if metric_2c is None:
        raise ValueError("DFTensor.metric not available")
    if ia_p is None:
        raise ValueError("DFTensor.iaP not available")
    if ij_p is None:
        raise ValueError("DFTensor.ijP not available")
    if cderi_3c is None:
        raise ValueError("DFTensor.cderi not available")

    # Ensure C-contiguous layout for efficient transfer
    metric_2c = np.ascontiguousarray(metric_2c, dtype=np.float64)
    cderi_3c = np.ascontiguousarray(cderi_3c, dtype=np.float64)
    ia_p = np.ascontiguousarray(ia_p, dtype=np.float64)
    ij_p = np.ascontiguousarray(ij_p, dtype=np.float64)

    # Get basis set name (handle different attribute names)
    basis = getattr(df_tensor, 'basis', None)
    if basis is None:
        # Try to get from molecule
        if hasattr(df_tensor, 'mol') and hasattr(df_tensor.mol, 'basis'):
            basis_obj = df_tensor.mol.basis
            # PySCF basis can be a string or dict
            if isinstance(basis_obj, str):
                basis = basis_obj
            elif isinstance(basis_obj, dict):
                # Get first basis set name from dict
                basis = next(iter(basis_obj.values())) if basis_obj else "unknown"
            else:
                basis = str(basis_obj)
        else:
            basis = "unknown"

    # Create Rust DFTensorsHDF5 object
    rust_df = DFTensorsHDF5(
        metric_2c=metric_2c,
        cderi_3c=cderi_3c,
        ia_p=ia_p,
        ij_p=ij_p,
        naux=df_tensor.naux,
        nao=df_tensor.nao,
        nocc=df_tensor.nocc,
        nvir=df_tensor.nvir,
        auxbasis=df_tensor.auxbasis,
        basis=basis,
    )

    # Write to HDF5
    rust_df.save_hdf5(filename)


def load_df_tensors_hdf5(filename: str):
    """
    Load DFTensor from HDF5 using Rust implementation

    Args:
        filename: Input HDF5 file path

    Returns:
        Dictionary containing loaded tensors and metadata
        Keys: 'metric_2c', 'cderi_3c', 'ia_p', 'ij_p', 'metadata'

    Raises:
        ImportError: If Rust HDF5 backend is not available
        IOError: If HDF5 read fails
    """
    if not HAS_RUST_HDF5:
        raise ImportError(
            "Rust HDF5 backend not available. "
            "Please rebuild with: maturin develop --release --features hdf5_support"
        )

    # Load from Rust
    rust_df = DFTensorsHDF5.load_hdf5(filename)

    # Extract data
    result = {
        'metric_2c': rust_df.get_metric_2c(),
        'cderi_3c': rust_df.get_cderi_3c(),
        'ia_p': rust_df.get_ia_p(),
        'ij_p': rust_df.get_ij_p(),
        'metadata': rust_df.get_metadata(),
    }

    return result


def convert_to_dftensor(loaded_data, mol):
    """
    Convert loaded HDF5 data back to a DFTensor object

    Args:
        loaded_data: Dictionary returned by load_df_tensors_hdf5()
        mol: PySCF molecule object

    Returns:
        DFTensor object
    """
    from quasix.df_tensor import DFTensor

    meta = loaded_data['metadata']

    # Create empty DFTensor
    df = DFTensor(mol, auxbasis=meta['auxbasis'])

    # Populate with loaded data
    df.metric = loaded_data['metric_2c']
    df.cderi = loaded_data['cderi_3c']
    df.iaP = loaded_data['ia_p']
    df.ijP = loaded_data['ij_p']

    # Set dimensions
    df.naux = meta['naux']
    df.nao = meta['nao']
    df.nocc = meta['nocc']
    df.nvir = meta['nvir']

    return df
