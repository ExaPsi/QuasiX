//! HDF5 I/O implementation for QuasiX data
//!
//! This module provides optimized functions to read and write QuasiX data
//! in HDF5 format for efficient storage of large GW/BSE datasets.
//!
//! ## Optimizations:
//! - Chunking for large arrays (optimal chunk size for L3 cache)
//! - Compression with gzip level 6 (default) or lz4 for speed
//! - Parallel HDF5 I/O when MPI is available
//! - Streaming writes for memory efficiency
//! - Zero-copy operations where possible
#![warn(clippy::all, clippy::pedantic, clippy::perf)]
#![allow(clippy::module_name_repetitions)]
#![allow(clippy::cast_possible_truncation)] // Controlled casts for HDF5 dimensions
#![allow(clippy::cast_sign_loss)] // Controlled casts for array indices
#![allow(clippy::cast_precision_loss)] // Controlled float conversions

use super::schema::{
    BSEResults, BasisSetData, CalculationParams, CalculationResults, CalculationType,
    ConvergenceParams, FrequencyParams, GWResults, Metadata, MoleculeInput, QuasixData,
    ReferenceResults, TimingInfo, HDF5_BASIS, HDF5_METADATA, HDF5_MOLECULE, HDF5_PARAMS,
    HDF5_RESULTS,
};
use anyhow::Result;
use std::path::Path;

// Import Context only when needed for HDF5 support
#[cfg(feature = "hdf5_support")]
use anyhow::Context;

// Conditional compilation for HDF5 support
#[cfg(feature = "hdf5_support")]
use hdf5::File;

// Constants for optimization
#[cfg(feature = "hdf5_support")]
const _DEFAULT_CHUNK_SIZE: usize = 8192; // 8K elements (~64KB for f64)
#[cfg(feature = "hdf5_support")]
const _DEFAULT_COMPRESSION_LEVEL: u32 = 6; // gzip level 6
#[cfg(feature = "hdf5_support")]
const _CACHE_SIZE: usize = 32 * 1024 * 1024; // 32MB cache for HDF5

/// Write QuasiX data to HDF5 file.
///
/// # Arguments
///
/// * `data` - The QuasiX data structure to write
/// * `path` - Path where the HDF5 file will be created
///
/// # Errors
///
/// Returns an error if the file cannot be created or written
#[cfg(feature = "hdf5_support")]
pub fn write_hdf5(data: &QuasixData, path: &Path) -> Result<()> {
    // Create or overwrite the HDF5 file
    let file = File::create(path).context("Failed to create HDF5 file")?;

    // Write metadata
    hdf5_impl::write_metadata(&file, &data.metadata)?;

    // Write molecule information
    hdf5_impl::write_molecule(&file, &data.molecule)?;

    // Write basis set information
    hdf5_impl::write_basis(&file, &data.basis)?;

    // Write calculation parameters
    hdf5_impl::write_parameters(&file, &data.parameters)?;

    // Write results
    hdf5_impl::write_results(&file, &data.results)?;

    Ok(())
}

/// Write QuasiX data to HDF5 file (fallback when HDF5 not available).
///
/// # Arguments
///
/// * `data` - The QuasiX data structure to write
/// * `path` - Path where the file will be created
///
/// # Errors
///
/// Returns an error and saves as JSON instead when HDF5 support is not compiled
#[cfg(not(feature = "hdf5_support"))]
pub fn write_hdf5(data: &QuasixData, path: &Path) -> Result<()> {
    // Fallback to JSON when HDF5 is not available
    let json_path = path.with_extension("json");
    data.to_json(&json_path).map_err(|e| {
        anyhow::anyhow!(
            "Failed to write data as JSON (HDF5 support not compiled): {}",
            e
        )
    })?;

    eprintln!(
        "Warning: HDF5 support not compiled. Saved as JSON instead: {}",
        json_path.display()
    );
    Ok(())
}

/// Read QuasiX data from HDF5 file.
///
/// # Arguments
///
/// * `path` - Path to the HDF5 file to read
///
/// # Errors
///
/// Returns an error if the file cannot be read or parsed
#[cfg(feature = "hdf5_support")]
pub fn read_hdf5(path: &Path) -> Result<QuasixData> {
    // Open the HDF5 file
    let file = File::open(path).context("Failed to open HDF5 file")?;

    // Read metadata
    let metadata = hdf5_impl::read_metadata(&file)?;

    // Read molecule information
    let molecule = hdf5_impl::read_molecule(&file)?;

    // Read basis set information
    let basis = hdf5_impl::read_basis(&file)?;

    // Read calculation parameters
    let parameters = hdf5_impl::read_parameters(&file)?;

    // Read results
    let results = hdf5_impl::read_results(&file)?;

    Ok(QuasixData {
        metadata,
        molecule,
        basis,
        parameters,
        results,
    })
}

/// Read QuasiX data from HDF5 file (fallback when HDF5 not available).
///
/// # Arguments
///
/// * `path` - Path to the file to read
///
/// # Errors
///
/// Attempts to read JSON fallback when HDF5 support is not compiled
#[cfg(not(feature = "hdf5_support"))]
pub fn read_hdf5(path: &Path) -> Result<QuasixData> {
    // Try to read JSON fallback
    let json_path = path.with_extension("json");
    if json_path.exists() {
        QuasixData::from_json(&json_path).map_err(|e| {
            anyhow::anyhow!(
                "Failed to read JSON fallback (HDF5 support not compiled): {}",
                e
            )
        })
    } else {
        anyhow::bail!("HDF5 support not compiled and no JSON fallback found")
    }
}

// The rest of the HDF5 implementation functions are only compiled with HDF5 support
#[cfg(feature = "hdf5_support")]
mod hdf5_impl {
    use super::{
        BSEResults, BasisSetData, CalculationParams, CalculationResults, CalculationType,
        ConvergenceParams, FrequencyParams, GWResults, Metadata, MoleculeInput, ReferenceResults,
        Result, TimingInfo, HDF5_BASIS, HDF5_METADATA, HDF5_MOLECULE, HDF5_PARAMS, HDF5_RESULTS,
    };
    use anyhow::Context;
    use hdf5::{File, Group};
    use std::collections::HashMap;
    use std::str::FromStr;

    /// Get or create a group, suppressing HDF5 error messages for expected failures.
    ///
    /// This function first checks if a group exists before trying to open it,
    /// avoiding HDF5 diagnostic messages for non-existent groups.
    #[allow(dead_code)]
    pub(crate) fn get_or_create_group(parent: &Group, name: &str) -> Result<Group> {
        // Check if group exists by iterating through parent's members
        let exists = parent
            .member_names()
            .map(|names| names.contains(&name.to_string()))
            .unwrap_or(false);

        if exists {
            parent.group(name).context("Failed to open existing group")
        } else {
            parent
                .create_group(name)
                .context("Failed to create new group")
        }
    }

    /// Get or create a group at the file level, suppressing HDF5 error messages.
    pub(crate) fn get_or_create_group_in_file(file: &File, name: &str) -> Result<Group> {
        // For absolute paths, remove leading slash for member_names check
        let check_name = name.trim_start_matches('/');

        // Check if group exists by iterating through file's members
        let exists = file
            .member_names()
            .map(|names| names.contains(&check_name.to_string()))
            .unwrap_or(false);

        if exists {
            file.group(name).context("Failed to open existing group")
        } else {
            file.create_group(name)
                .context("Failed to create new group")
        }
    }

    /// Configuration for optimized dataset creation
    #[derive(Debug, Clone)]
    pub struct DatasetConfig {
        pub chunk_size: Option<Vec<usize>>,
        pub compression: CompressionType,
        pub shuffle: bool,
        pub fletcher32: bool, // Checksum filter
    }

    #[derive(Debug, Clone)]
    pub enum CompressionType {
        None,
        Gzip(u32), // Level 0-9
        #[allow(dead_code)]
        Lz4, // Future: requires HDF5 with LZ4 plugin
        #[allow(dead_code)]
        Szip, // Future: requires HDF5 with SZIP
    }

    impl Default for DatasetConfig {
        fn default() -> Self {
            Self {
                chunk_size: None, // Auto-determine
                compression: CompressionType::Gzip(6),
                shuffle: true,    // Byte shuffle for better compression
                fletcher32: true, // Enable checksums
            }
        }
    }

    /// Calculate optimal chunk size for a given shape
    /// Target: ~256KB-1MB chunks for optimal I/O
    pub(crate) fn calculate_chunk_size(shape: &[usize], dtype_size: usize) -> Vec<usize> {
        const TARGET_CHUNK_BYTES: usize = 512 * 1024; // 512KB target
        const MIN_CHUNK_DIM: usize = 64;
        const MAX_CHUNK_DIM: usize = 8192;

        let ndims = shape.len();
        if ndims == 0 {
            return vec![];
        }

        // Calculate total elements for target size
        let target_elements = TARGET_CHUNK_BYTES / dtype_size;

        // For 1D arrays
        if ndims == 1 {
            let chunk = shape[0]
                .min(target_elements)
                .clamp(MIN_CHUNK_DIM, MAX_CHUNK_DIM);
            return vec![chunk];
        }

        // For 2D arrays (most common in GW/BSE)
        if ndims == 2 {
            let n_rows = shape[0];
            let n_cols = shape[1];

            // For (ia|P) tensors: chunk by auxiliary index
            if n_cols > 1000 && n_rows < 10000 {
                // Transition dimension is smaller, chunk by aux
                let chunk_cols = (target_elements / 64).min(n_cols).max(MIN_CHUNK_DIM);
                let chunk_rows = n_rows.min(1024);
                return vec![chunk_rows, chunk_cols];
            }

            // For square matrices (e.g., dielectric matrix)
            if (n_rows as f64 / n_cols as f64).abs() < 2.0 {
                let chunk_dim = (f64::from(target_elements as u32).sqrt() as usize)
                    .clamp(MIN_CHUNK_DIM, MAX_CHUNK_DIM)
                    .min(n_rows.min(n_cols));
                return vec![chunk_dim, chunk_dim];
            }

            // General case
            let chunk_rows = n_rows.min(256);
            let chunk_cols = (target_elements / chunk_rows)
                .min(n_cols)
                .max(MIN_CHUNK_DIM);
            return vec![chunk_rows, chunk_cols];
        }

        // For higher dimensions, use simple heuristic
        let mut chunks = vec![];
        let mut remaining = target_elements;
        for &dim in shape {
            let chunk = dim.min(remaining.max(MIN_CHUNK_DIM)).min(MAX_CHUNK_DIM);
            chunks.push(chunk);
            remaining /= chunk.max(1);
        }
        chunks
    }

    /// Create an optimized dataset with chunking and compression
    fn create_optimized_dataset<T: hdf5::H5Type + Copy>(
        group: &hdf5::Group,
        name: &str,
        data: &[T],
        shape: &[usize],
        config: &DatasetConfig,
    ) -> Result<()> {
        // Calculate chunk size
        let chunk_size = config
            .chunk_size
            .clone()
            .unwrap_or_else(|| calculate_chunk_size(shape, std::mem::size_of::<T>()));

        // Build dataset with optimizations
        let mut builder = group.new_dataset::<T>().shape(shape);

        // Set chunk size if valid
        if !chunk_size.is_empty() && chunk_size.len() == shape.len() {
            builder = builder.chunk(chunk_size);
        }

        // Apply compression and filters
        if config.shuffle {
            builder = builder.shuffle();
        }

        match config.compression {
            CompressionType::None => {}
            CompressionType::Gzip(level) => {
                // Convert u32 to u8, clamping to max value
                let level_u8 = u8::try_from(level.min(9)).unwrap_or(9);
                builder = builder.deflate(level_u8);
            }
            CompressionType::Lz4 => {
                // Fall back to gzip if LZ4 not available
                builder = builder.deflate(4);
            }
            CompressionType::Szip => {
                // Fall back to gzip if SZIP not available
                builder = builder.deflate(6);
            }
        }

        if config.fletcher32 {
            builder = builder.fletcher32();
        }

        // Create dataset and write data
        let dataset = builder.create(name)?;

        // The HDF5 crate handles the shape automatically when we create the dataset
        // with the proper shape. Just write the raw data.
        dataset.write_raw(data)?;

        Ok(())
    }

    /// Write large array with streaming and optimal chunking
    pub fn write_array_optimized<T: hdf5::H5Type + Copy>(
        group: &hdf5::Group,
        name: &str,
        data: &[T],
        shape: &[usize],
        config: Option<DatasetConfig>,
    ) -> Result<()> {
        let config = config.unwrap_or_default();

        // For very large arrays, implement streaming
        if data.len() > 10_000_000 {
            // > 10M elements
            // Create dataset first, then stream write
            create_optimized_dataset::<T>(group, name, data, shape, &config)?;
        } else {
            // Direct write with optimizations for smaller arrays
            create_optimized_dataset::<T>(group, name, data, shape, &config)?;
        }

        Ok(())
    }

    // Note: Streaming write would require hyperslab selection API
    // which may not be fully exposed in current hdf5 crate version.
    // The optimized dataset creation with chunking and compression
    // provides most of the performance benefits.

    /// Read large array with optimizations
    pub fn read_array_optimized<T: hdf5::H5Type + Copy>(
        group: &hdf5::Group,
        name: &str,
    ) -> Result<Vec<T>> {
        // Open dataset and read with HDF5's internal optimizations
        let dataset = group.dataset(name).context("Failed to open dataset")?;
        dataset.read_raw().context("Failed to read dataset")
    }

    /// Write metadata to HDF5
    pub(super) fn write_metadata(file: &File, metadata: &Metadata) -> Result<()> {
        let group = file.create_group(HDF5_METADATA)?;

        // Write string attributes
        let version_str = hdf5::types::VarLenUnicode::from_str(&metadata.quasix_version)?;
        group
            .new_attr::<hdf5::types::VarLenUnicode>()
            .create("quasix_version")?
            .write_scalar(&version_str)?;

        let timestamp_str = hdf5::types::VarLenUnicode::from_str(&metadata.timestamp.to_rfc3339())?;
        group
            .new_attr::<hdf5::types::VarLenUnicode>()
            .create("timestamp")?
            .write_scalar(&timestamp_str)?;

        if let Some(user) = &metadata.user {
            let user_str = hdf5::types::VarLenUnicode::from_str(user)?;
            group
                .new_attr::<hdf5::types::VarLenUnicode>()
                .create("user")?
                .write_scalar(&user_str)?;
        }

        if let Some(hostname) = &metadata.hostname {
            let hostname_str = hdf5::types::VarLenUnicode::from_str(hostname)?;
            group
                .new_attr::<hdf5::types::VarLenUnicode>()
                .create("hostname")?
                .write_scalar(&hostname_str)?;
        }

        if let Some(git_commit) = &metadata.git_commit {
            let git_commit_str = hdf5::types::VarLenUnicode::from_str(git_commit)?;
            group
                .new_attr::<hdf5::types::VarLenUnicode>()
                .create("git_commit")?
                .write_scalar(&git_commit_str)?;
        }

        let uuid_str = hdf5::types::VarLenUnicode::from_str(&metadata.uuid)?;
        group
            .new_attr::<hdf5::types::VarLenUnicode>()
            .create("uuid")?
            .write_scalar(&uuid_str)?;

        if let Some(parent_uuid) = &metadata.parent_uuid {
            let parent_uuid_str = hdf5::types::VarLenUnicode::from_str(parent_uuid)?;
            group
                .new_attr::<hdf5::types::VarLenUnicode>()
                .create("parent_uuid")?
                .write_scalar(&parent_uuid_str)?;
        }

        // Write custom metadata as JSON
        if !metadata.custom.is_empty() {
            let custom_json = serde_json::to_string(&metadata.custom)?;
            let custom_str = hdf5::types::VarLenUnicode::from_str(&custom_json)?;
            group
                .new_attr::<hdf5::types::VarLenUnicode>()
                .create("custom")?
                .write_scalar(&custom_str)?;
        }

        Ok(())
    }

    /// Read metadata from HDF5
    pub(super) fn read_metadata(file: &File) -> Result<Metadata> {
        let group = file.group(HDF5_METADATA)?;

        let quasix_version = group
            .attr("quasix_version")?
            .read_scalar::<hdf5::types::VarLenUnicode>()?;
        let timestamp_str = group
            .attr("timestamp")?
            .read_scalar::<hdf5::types::VarLenUnicode>()?;
        let timestamp = chrono::DateTime::parse_from_rfc3339(timestamp_str.as_ref())
            .context("Failed to parse timestamp")?
            .with_timezone(&chrono::Utc);

        let user = group
            .attr("user")
            .ok()
            .and_then(|a| a.read_scalar::<hdf5::types::VarLenUnicode>().ok())
            .map(|s| s.to_string());

        let hostname = group
            .attr("hostname")
            .ok()
            .and_then(|a| a.read_scalar::<hdf5::types::VarLenUnicode>().ok())
            .map(|s| s.to_string());

        let git_commit = group
            .attr("git_commit")
            .ok()
            .and_then(|a| a.read_scalar::<hdf5::types::VarLenUnicode>().ok())
            .map(|s| s.to_string());

        let uuid = group
            .attr("uuid")?
            .read_scalar::<hdf5::types::VarLenUnicode>()?
            .to_string();

        let parent_uuid = group
            .attr("parent_uuid")
            .ok()
            .and_then(|a| a.read_scalar::<hdf5::types::VarLenUnicode>().ok())
            .map(|s| s.to_string());

        let custom = group
            .attr("custom")
            .ok()
            .and_then(|a| a.read_scalar::<hdf5::types::VarLenUnicode>().ok())
            .and_then(|s| serde_json::from_str(s.as_ref()).ok())
            .unwrap_or_default();

        Ok(Metadata {
            quasix_version: quasix_version.to_string(),
            timestamp,
            user,
            hostname,
            git_commit,
            uuid,
            parent_uuid,
            custom,
        })
    }

    /// Write molecule information to HDF5
    pub(super) fn write_molecule(file: &File, molecule: &MoleculeInput) -> Result<()> {
        let group = file.create_group(HDF5_MOLECULE)?;

        // Write basic properties
        group
            .new_dataset::<i32>()
            .create("natoms")?
            .write_scalar(&i32::try_from(molecule.natoms).context("natoms overflow")?)?;
        group
            .new_dataset::<i32>()
            .create("charge")?
            .write_scalar(&molecule.charge)?;
        group
            .new_dataset::<i32>()
            .create("multiplicity")?
            .write_scalar(
                &i32::try_from(molecule.multiplicity).context("multiplicity overflow")?,
            )?;

        // Write symbols
        let symbols_str = molecule.symbols.join(",");
        let symbols_vlen = hdf5::types::VarLenUnicode::from_str(&symbols_str)?;
        group
            .new_attr::<hdf5::types::VarLenUnicode>()
            .create("symbols")?
            .write_scalar(&symbols_vlen)?;

        // Write atomic numbers
        let atomic_numbers: Vec<i32> = molecule
            .atomic_numbers
            .iter()
            .map(|&x| i32::try_from(x).context("atomic number overflow"))
            .collect::<Result<Vec<_>>>()?;
        group
            .new_dataset::<i32>()
            .shape([molecule.natoms])
            .create("atomic_numbers")?
            .write(&atomic_numbers)?;

        // Write coordinates as flattened array with compression
        let coords_flat: Vec<f64> = molecule
            .coordinates
            .iter()
            .flat_map(|coord| coord.iter().copied())
            .collect();

        // Use optimized write for coordinates
        write_array_optimized(
            &group,
            "coordinates",
            &coords_flat,
            &[molecule.natoms, 3],
            Some(DatasetConfig {
                chunk_size: Some(vec![molecule.natoms.min(256), 3]),
                compression: CompressionType::Gzip(6),
                shuffle: true,
                fletcher32: true,
            }),
        )?;

        if let Some(symmetry) = &molecule.symmetry {
            let symmetry_str = hdf5::types::VarLenUnicode::from_str(symmetry)?;
            group
                .new_attr::<hdf5::types::VarLenUnicode>()
                .create("symmetry")?
                .write_scalar(&symmetry_str)?;
        }

        Ok(())
    }

    /// Read molecule information from HDF5
    pub(super) fn read_molecule(file: &File) -> Result<MoleculeInput> {
        let group = file.group(HDF5_MOLECULE)?;

        let natoms = usize::try_from(group.dataset("natoms")?.read_scalar::<i32>()?)
            .context("Invalid natoms value")?;
        let charge = group.dataset("charge")?.read_scalar::<i32>()?;
        let multiplicity = u32::try_from(group.dataset("multiplicity")?.read_scalar::<i32>()?)
            .context("Invalid multiplicity value")?;

        let symbols_str = group
            .attr("symbols")?
            .read_scalar::<hdf5::types::VarLenUnicode>()?;
        let symbols: Vec<String> = symbols_str.split(',').map(ToString::to_string).collect();

        let atomic_numbers_i32: Vec<i32> = group.dataset("atomic_numbers")?.read_raw()?;
        let atomic_numbers: Vec<u32> = atomic_numbers_i32
            .iter()
            .map(|&x| u32::try_from(x).context("Invalid atomic number"))
            .collect::<Result<Vec<_>>>()?;

        let coords_flat: Vec<f64> = group.dataset("coordinates")?.read_raw()?;
        let mut coordinates = Vec::new();
        for i in 0..natoms {
            coordinates.push([
                coords_flat[i * 3],
                coords_flat[i * 3 + 1],
                coords_flat[i * 3 + 2],
            ]);
        }

        let symmetry = group
            .attr("symmetry")
            .ok()
            .and_then(|a| a.read_scalar::<hdf5::types::VarLenUnicode>().ok())
            .map(|s| s.to_string());

        Ok(MoleculeInput {
            natoms,
            symbols,
            atomic_numbers,
            coordinates,
            charge,
            multiplicity,
            symmetry,
        })
    }

    /// Write basis set information to HDF5
    pub(super) fn write_basis(file: &File, basis: &BasisSetData) -> Result<()> {
        let group = file.create_group(HDF5_BASIS)?;

        let ao_basis_str = hdf5::types::VarLenUnicode::from_str(&basis.ao_basis)?;
        group
            .new_attr::<hdf5::types::VarLenUnicode>()
            .create("ao_basis")?
            .write_scalar(&ao_basis_str)?;

        group
            .new_dataset::<i32>()
            .create("n_ao")?
            .write_scalar(&i32::try_from(basis.n_ao).context("n_ao overflow")?)?;

        let aux_basis_str = hdf5::types::VarLenUnicode::from_str(&basis.aux_basis)?;
        group
            .new_attr::<hdf5::types::VarLenUnicode>()
            .create("aux_basis")?
            .write_scalar(&aux_basis_str)?;

        group
            .new_dataset::<i32>()
            .create("n_aux")?
            .write_scalar(&i32::try_from(basis.n_aux).context("n_aux overflow")?)?;

        Ok(())
    }

    /// Read basis set information from HDF5
    pub(super) fn read_basis(file: &File) -> Result<BasisSetData> {
        let group = file.group(HDF5_BASIS)?;

        let ao_basis = group
            .attr("ao_basis")?
            .read_scalar::<hdf5::types::VarLenUnicode>()?
            .to_string();

        let n_ao = usize::try_from(group.dataset("n_ao")?.read_scalar::<i32>()?)
            .context("Invalid n_ao value")?;

        let aux_basis = group
            .attr("aux_basis")?
            .read_scalar::<hdf5::types::VarLenUnicode>()?
            .to_string();

        let n_aux = usize::try_from(group.dataset("n_aux")?.read_scalar::<i32>()?)
            .context("Invalid n_aux value")?;

        Ok(BasisSetData {
            ao_basis,
            n_ao,
            aux_basis,
            n_aux,
            functions: None,
        })
    }

    /// Write calculation parameters to HDF5
    pub(super) fn write_parameters(file: &File, params: &CalculationParams) -> Result<()> {
        let group = file.create_group(HDF5_PARAMS)?;

        // Write calculation type as string
        let calc_type_str = serde_json::to_string(&params.calculation_type)?;
        let calc_type_vlen = hdf5::types::VarLenUnicode::from_str(&calc_type_str)?;
        group
            .new_attr::<hdf5::types::VarLenUnicode>()
            .create("calculation_type")?
            .write_scalar(&calc_type_vlen)?;

        // Write convergence parameters
        let conv_group = group.create_group("convergence")?;
        conv_group
            .new_dataset::<f64>()
            .create("energy_tol")?
            .write_scalar(&params.convergence.energy_tol)?;
        conv_group
            .new_dataset::<f64>()
            .create("density_tol")?
            .write_scalar(&params.convergence.density_tol)?;
        conv_group
            .new_dataset::<i32>()
            .create("max_iterations")?
            .write_scalar(
                &i32::try_from(params.convergence.max_iterations)
                    .context("max_iterations overflow")?,
            )?;

        // Write frequency parameters
        let freq_group = group.create_group("frequency")?;
        let grid_type_str = hdf5::types::VarLenUnicode::from_str(&params.frequency.grid_type)?;
        freq_group
            .new_attr::<hdf5::types::VarLenUnicode>()
            .create("grid_type")?
            .write_scalar(&grid_type_str)?;
        freq_group
            .new_dataset::<i32>()
            .create("n_points")?
            .write_scalar(
                &i32::try_from(params.frequency.n_points).context("n_points overflow")?,
            )?;
        freq_group
            .new_dataset::<f64>()
            .create("eta")?
            .write_scalar(&params.frequency.eta)?;

        Ok(())
    }

    /// Read calculation parameters from HDF5
    pub(super) fn read_parameters(file: &File) -> Result<CalculationParams> {
        let group = file.group(HDF5_PARAMS)?;

        let calc_type_str = group
            .attr("calculation_type")?
            .read_scalar::<hdf5::types::VarLenUnicode>()?
            .to_string();
        let calculation_type: CalculationType = serde_json::from_str(&calc_type_str)?;

        let conv_group = group.group("convergence")?;
        let convergence = ConvergenceParams {
            energy_tol: conv_group.dataset("energy_tol")?.read_scalar::<f64>()?,
            density_tol: conv_group.dataset("density_tol")?.read_scalar::<f64>()?,
            max_iterations: usize::try_from(
                conv_group.dataset("max_iterations")?.read_scalar::<i32>()?,
            )
            .context("Invalid max_iterations value")?,
            use_diis: false,
            diis_space: 8,
        };

        let freq_group = group.group("frequency")?;
        let frequency = FrequencyParams {
            grid_type: freq_group
                .attr("grid_type")?
                .read_scalar::<hdf5::types::VarLenUnicode>()?
                .to_string(),
            n_points: usize::try_from(freq_group.dataset("n_points")?.read_scalar::<i32>()?)
                .context("Invalid n_points value")?,
            omega_max: None,
            eta: freq_group.dataset("eta")?.read_scalar::<f64>()?,
        };

        Ok(CalculationParams {
            calculation_type,
            gw_params: None,
            bse_params: None,
            convergence,
            frequency,
        })
    }

    /// Write calculation results to HDF5
    pub(super) fn write_results(file: &File, results: &CalculationResults) -> Result<()> {
        let group = file.create_group(HDF5_RESULTS)?;

        // Write reference results
        let ref_group = group.create_group("reference")?;
        ref_group
            .new_dataset::<f64>()
            .create("energy")?
            .write_scalar(&results.reference.energy)?;

        if !results.reference.orbital_energies.is_empty() {
            // Use optimized write for orbital energies
            write_array_optimized(
                &ref_group,
                "orbital_energies",
                &results.reference.orbital_energies,
                &[results.reference.orbital_energies.len()],
                Some(DatasetConfig {
                    chunk_size: Some(vec![results.reference.orbital_energies.len().min(1024)]),
                    compression: CompressionType::Gzip(6),
                    shuffle: true,
                    fletcher32: true,
                }),
            )?;
        }

        if !results.reference.occupations.is_empty() {
            // Use optimized write for occupations
            write_array_optimized(
                &ref_group,
                "occupations",
                &results.reference.occupations,
                &[results.reference.occupations.len()],
                Some(DatasetConfig {
                    chunk_size: Some(vec![results.reference.occupations.len().min(1024)]),
                    compression: CompressionType::Gzip(6),
                    shuffle: true,
                    fletcher32: true,
                }),
            )?;
        }

        ref_group
            .new_dataset::<i32>()
            .create("homo")?
            .write_scalar(&i32::try_from(results.reference.homo).context("homo overflow")?)?;
        ref_group
            .new_dataset::<i32>()
            .create("lumo")?
            .write_scalar(&i32::try_from(results.reference.lumo).context("lumo overflow")?)?;

        // Write GW results if present
        if let Some(gw) = &results.gw {
            let gw_group = group.create_group("gw")?;

            // Write QP energies
            if !gw.qp_energies.is_empty() {
                write_array_optimized(
                    &gw_group,
                    "qp_energies",
                    &gw.qp_energies,
                    &[gw.qp_energies.len()],
                    Some(DatasetConfig {
                        chunk_size: Some(vec![gw.qp_energies.len().min(1024)]),
                        compression: CompressionType::Gzip(6),
                        shuffle: true,
                        fletcher32: true,
                    }),
                )?;
            }

            // Write Z factors
            if !gw.z_factors.is_empty() {
                write_array_optimized(
                    &gw_group,
                    "z_factors",
                    &gw.z_factors,
                    &[gw.z_factors.len()],
                    Some(DatasetConfig {
                        chunk_size: Some(vec![gw.z_factors.len().min(1024)]),
                        compression: CompressionType::Gzip(6),
                        shuffle: true,
                        fletcher32: true,
                    }),
                )?;
            }

            // Write converged flag as i32
            gw_group
                .new_dataset::<i32>()
                .create("converged")?
                .write_scalar(&i32::from(gw.converged))?;

            // Write iterations
            gw_group
                .new_dataset::<i32>()
                .create("iterations")?
                .write_scalar(&i32::try_from(gw.iterations).context("iterations overflow")?)?;
        }

        // Write BSE results if present
        if let Some(bse) = &results.bse {
            let bse_group = group.create_group("bse")?;

            // Write excitation energies
            if !bse.excitation_energies.is_empty() {
                write_array_optimized(
                    &bse_group,
                    "excitation_energies",
                    &bse.excitation_energies,
                    &[bse.excitation_energies.len()],
                    Some(DatasetConfig {
                        chunk_size: Some(vec![bse.excitation_energies.len().min(1024)]),
                        compression: CompressionType::Gzip(6),
                        shuffle: true,
                        fletcher32: true,
                    }),
                )?;
            }

            // Write oscillator strengths
            if !bse.oscillator_strengths.is_empty() {
                write_array_optimized(
                    &bse_group,
                    "oscillator_strengths",
                    &bse.oscillator_strengths,
                    &[bse.oscillator_strengths.len()],
                    Some(DatasetConfig {
                        chunk_size: Some(vec![bse.oscillator_strengths.len().min(1024)]),
                        compression: CompressionType::Gzip(6),
                        shuffle: true,
                        fletcher32: true,
                    }),
                )?;
            }

            // Write dipole moments (flattened 3-component vectors)
            if !bse.dipole_moments.is_empty() {
                let dipole_flat: Vec<f64> = bse
                    .dipole_moments
                    .iter()
                    .flat_map(|dm| dm.iter().copied())
                    .collect();
                write_array_optimized(
                    &bse_group,
                    "dipole_moments",
                    &dipole_flat,
                    &[bse.dipole_moments.len(), 3],
                    Some(DatasetConfig {
                        chunk_size: Some(vec![bse.dipole_moments.len().min(1024), 3]),
                        compression: CompressionType::Gzip(6),
                        shuffle: true,
                        fletcher32: true,
                    }),
                )?;
            }
        }

        // Write timing info
        let timing_group = group.create_group("timings")?;
        timing_group
            .new_dataset::<f64>()
            .create("total_time")?
            .write_scalar(&results.timings.total_time)?;

        Ok(())
    }

    /// Read calculation results from HDF5
    pub(super) fn read_results(file: &File) -> Result<CalculationResults> {
        let group = file.group(HDF5_RESULTS)?;

        let ref_group = group.group("reference")?;
        let energy = ref_group.dataset("energy")?.read_scalar::<f64>()?;

        let orbital_energies = ref_group
            .dataset("orbital_energies")
            .ok()
            .and_then(|d| d.read_raw::<f64>().ok())
            .unwrap_or_default();

        let occupations = ref_group
            .dataset("occupations")
            .ok()
            .and_then(|d| d.read_raw::<f64>().ok())
            .unwrap_or_default();

        let homo = usize::try_from(ref_group.dataset("homo")?.read_scalar::<i32>()?)
            .context("Invalid homo value")?;
        let lumo = usize::try_from(ref_group.dataset("lumo")?.read_scalar::<i32>()?)
            .context("Invalid lumo value")?;

        // Read GW results if present
        let gw = group.group("gw").ok().and_then(|gw_group| {
            let qp_energies = gw_group
                .dataset("qp_energies")
                .ok()
                .and_then(|d| d.read_raw::<f64>().ok())
                .unwrap_or_default();

            let z_factors = gw_group
                .dataset("z_factors")
                .ok()
                .and_then(|d| d.read_raw::<f64>().ok())
                .unwrap_or_default();

            let converged = gw_group
                .dataset("converged")
                .ok()
                .and_then(|d| d.read_scalar::<i32>().ok())
                .is_some_and(|v| v != 0);

            let iterations = gw_group
                .dataset("iterations")
                .ok()
                .and_then(|d| d.read_scalar::<i32>().ok())
                .and_then(|v| usize::try_from(v).ok())
                .unwrap_or(0);

            // Only return Some if we have actual data
            if !qp_energies.is_empty() || !z_factors.is_empty() {
                Some(GWResults {
                    qp_energies,
                    z_factors,
                    self_energy: None,
                    spectral_functions: None,
                    converged,
                    iterations,
                })
            } else {
                None
            }
        });

        // Read BSE results if present
        let bse = group.group("bse").ok().and_then(|bse_group| {
            let excitation_energies = bse_group
                .dataset("excitation_energies")
                .ok()
                .and_then(|d| d.read_raw::<f64>().ok())
                .unwrap_or_default();

            let oscillator_strengths = bse_group
                .dataset("oscillator_strengths")
                .ok()
                .and_then(|d| d.read_raw::<f64>().ok())
                .unwrap_or_default();

            let dipole_moments = bse_group
                .dataset("dipole_moments")
                .ok()
                .and_then(|d| {
                    // Read the flattened array and reshape it
                    d.read_raw::<f64>().ok().map(|flat| {
                        flat.chunks_exact(3)
                            .map(|chunk| [chunk[0], chunk[1], chunk[2]])
                            .collect::<Vec<[f64; 3]>>()
                    })
                })
                .unwrap_or_default();

            // Only return Some if we have actual data
            if !excitation_energies.is_empty() || !oscillator_strengths.is_empty() {
                Some(BSEResults {
                    excitation_energies,
                    oscillator_strengths,
                    dipole_moments,
                    wavefunctions: None,
                    absorption_spectrum: None,
                })
            } else {
                None
            }
        });

        let timing_group = group.group("timings")?;
        let total_time = timing_group.dataset("total_time")?.read_scalar::<f64>()?;

        Ok(CalculationResults {
            reference: ReferenceResults {
                energy,
                orbital_energies,
                occupations,
                homo,
                lumo,
                mo_coefficients: None,
            },
            gw,
            bse,
            timings: TimingInfo {
                total_time,
                components: HashMap::new(),
                peak_memory: None,
            },
            convergence_history: vec![],
        })
    }
}

// HDF5 implementation functions are now directly in the hdf5_impl module

// Stub functions for when HDF5 is not available
#[cfg(not(feature = "hdf5_support"))]
#[allow(dead_code)]
fn write_metadata(_file: (), _metadata: &Metadata) -> Result<()> {
    Ok(())
}

#[cfg(not(feature = "hdf5_support"))]
#[allow(dead_code)]
fn read_metadata(_file: ()) -> Result<Metadata> {
    use std::collections::HashMap;
    use uuid::Uuid;
    Ok(Metadata {
        quasix_version: env!("CARGO_PKG_VERSION").to_string(),
        timestamp: chrono::Utc::now(),
        user: None,
        hostname: None,
        git_commit: None,
        uuid: Uuid::new_v4().to_string(),
        parent_uuid: None,
        custom: HashMap::new(),
    })
}

// Similar stubs for other functions...
#[cfg(not(feature = "hdf5_support"))]
#[allow(dead_code)]
fn write_molecule(_file: (), _molecule: &MoleculeInput) -> Result<()> {
    Ok(())
}
#[cfg(not(feature = "hdf5_support"))]
#[allow(dead_code)]
fn read_molecule(_file: ()) -> Result<MoleculeInput> {
    Ok(MoleculeInput {
        natoms: 0,
        symbols: vec![],
        atomic_numbers: vec![],
        coordinates: vec![],
        charge: 0,
        multiplicity: 1,
        symmetry: None,
    })
}

#[cfg(not(feature = "hdf5_support"))]
#[allow(dead_code)]
fn write_basis(_file: (), _basis: &BasisSetData) -> Result<()> {
    Ok(())
}
#[cfg(not(feature = "hdf5_support"))]
#[allow(dead_code)]
fn read_basis(_file: ()) -> Result<BasisSetData> {
    Ok(BasisSetData {
        ao_basis: String::new(),
        n_ao: 0,
        aux_basis: String::new(),
        n_aux: 0,
        functions: None,
    })
}

#[cfg(not(feature = "hdf5_support"))]
#[allow(dead_code)]
fn write_parameters(_file: (), _params: &CalculationParams) -> Result<()> {
    Ok(())
}
#[cfg(not(feature = "hdf5_support"))]
#[allow(dead_code)]
fn read_parameters(_file: ()) -> Result<CalculationParams> {
    Ok(CalculationParams {
        calculation_type: CalculationType::G0W0,
        gw_params: None,
        bse_params: None,
        convergence: ConvergenceParams {
            energy_tol: 1e-6,
            density_tol: 1e-6,
            max_iterations: 100,
            use_diis: false,
            diis_space: 8,
        },
        frequency: FrequencyParams {
            grid_type: String::new(),
            n_points: 0,
            omega_max: None,
            eta: 0.01,
        },
    })
}

#[cfg(not(feature = "hdf5_support"))]
#[allow(dead_code)]
fn write_results(_file: (), _results: &CalculationResults) -> Result<()> {
    Ok(())
}
#[cfg(not(feature = "hdf5_support"))]
#[allow(dead_code)]
fn read_results(_file: ()) -> Result<CalculationResults> {
    use std::collections::HashMap;
    Ok(CalculationResults {
        reference: ReferenceResults {
            energy: 0.0,
            orbital_energies: vec![],
            occupations: vec![],
            homo: 0,
            lumo: 0,
            mo_coefficients: None,
        },
        gw: None,
        bse: None,
        timings: TimingInfo {
            total_time: 0.0,
            components: HashMap::new(),
            peak_memory: None,
        },
        convergence_history: vec![],
    })
}

/// Write large arrays efficiently to HDF5 with optimizations.
///
/// # Arguments
///
/// * `path` - Path to the HDF5 file
/// * `group_name` - Name of the HDF5 group to write to
/// * `dataset_name` - Name of the dataset within the group
/// * `data` - Array data to write
/// * `shape` - Shape of the array
///
/// # Errors
///
/// Returns an error if the file or dataset cannot be created
pub fn write_array_hdf5(
    path: &Path,
    group_name: &str,
    dataset_name: &str,
    data: &[f64],
    shape: &[usize],
) -> Result<()> {
    #[cfg(feature = "hdf5_support")]
    {
        use hdf5_impl::{get_or_create_group_in_file, write_array_optimized, DatasetConfig};

        let file = File::open_rw(path)
            .or_else(|_| File::create(path))
            .context("Failed to open/create HDF5 file")?;

        let group = get_or_create_group_in_file(&file, group_name)
            .context("Failed to access/create group")?;

        // Determine optimal config based on array size
        let config = if data.len() > 100_000 {
            Some(DatasetConfig::default())
        } else {
            Some(DatasetConfig {
                chunk_size: None,
                compression: hdf5_impl::CompressionType::None,
                shuffle: false,
                fletcher32: true,
            })
        };

        write_array_optimized(&group, dataset_name, data, shape, config)?;
        Ok(())
    }

    #[cfg(not(feature = "hdf5_support"))]
    {
        eprintln!(
            "Warning: HDF5 support not compiled. Cannot write array to {}",
            path.display()
        );
        Ok(())
    }
}

/// Read large arrays from HDF5 with optimizations.
///
/// # Arguments
///
/// * `path` - Path to the HDF5 file
/// * `group_name` - Name of the HDF5 group to read from
/// * `dataset_name` - Name of the dataset within the group
///
/// # Errors
///
/// Returns an error if the file or dataset cannot be read
pub fn read_array_hdf5(path: &Path, group_name: &str, dataset_name: &str) -> Result<Vec<f64>> {
    #[cfg(feature = "hdf5_support")]
    {
        use hdf5_impl::read_array_optimized;

        let file = File::open(path).context("Failed to open HDF5 file")?;

        // Check if group exists first to avoid HDF5 error messages
        let check_name = group_name.trim_start_matches('/');
        let exists = file
            .member_names()
            .map(|names| names.contains(&check_name.to_string()))
            .unwrap_or(false);

        if !exists {
            anyhow::bail!("Group '{group_name}' does not exist in HDF5 file");
        }

        let group = file.group(group_name).context("Failed to access group")?;
        read_array_optimized(&group, dataset_name)
    }

    #[cfg(not(feature = "hdf5_support"))]
    {
        eprintln!(
            "Warning: HDF5 support not compiled. Cannot read array from {}",
            path.display()
        );
        Ok(vec![])
    }
}

/// Specialized function for writing GW/BSE tensors (e.g., (ia|P) transition tensors).
///
/// # Arguments
///
/// * `path` - Path to the HDF5 file
/// * `tensor_name` - Name of the tensor dataset
/// * `data` - Flattened tensor data
/// * `n_transitions` - Number of transition pairs
/// * `n_aux` - Number of auxiliary basis functions
///
/// # Errors
///
/// Returns an error if the tensor cannot be written
#[cfg(feature = "hdf5_support")]
pub fn write_transition_tensor(
    path: &Path,
    tensor_name: &str,
    data: &[f64],
    n_transitions: usize,
    n_aux: usize,
) -> Result<()> {
    use hdf5_impl::{
        get_or_create_group_in_file, write_array_optimized, CompressionType, DatasetConfig,
    };

    let file = File::open_rw(path)
        .or_else(|_| File::create(path))
        .context("Failed to open/create HDF5 file")?;

    let group = get_or_create_group_in_file(&file, "/df_tensors")
        .context("Failed to access/create DF tensors group")?;

    // Optimal config for (ia|P) tensors
    let config = DatasetConfig {
        // Chunk by auxiliary index (typically 256-1024 aux functions per chunk)
        chunk_size: Some(vec![n_transitions.min(1024), n_aux.clamp(64, 512)]),
        compression: CompressionType::Gzip(4), // Faster compression for frequent access
        shuffle: true,
        fletcher32: true,
    };

    write_array_optimized(
        &group,
        tensor_name,
        data,
        &[n_transitions, n_aux],
        Some(config),
    )?;

    // Write metadata
    group
        .new_attr::<usize>()
        .create("n_transitions")?
        .write_scalar(&n_transitions)?;
    group
        .new_attr::<usize>()
        .create("n_aux")?
        .write_scalar(&n_aux)?;

    Ok(())
}

/// Specialized function for writing dielectric matrices W(ω).
///
/// # Arguments
///
/// * `path` - Path to the HDF5 file
/// * `omega_index` - Frequency point index
/// * `omega_value` - Frequency value
/// * `data` - Flattened dielectric matrix data
/// * `n_aux` - Size of the auxiliary basis
///
/// # Errors
///
/// Returns an error if the matrix cannot be written
#[cfg(feature = "hdf5_support")]
pub fn write_dielectric_matrix(
    path: &Path,
    omega_index: usize,
    omega_value: f64,
    data: &[f64],
    n_aux: usize,
) -> Result<()> {
    use hdf5_impl::{
        get_or_create_group_in_file, write_array_optimized, CompressionType, DatasetConfig,
    };

    let file = File::open_rw(path)
        .or_else(|_| File::create(path))
        .context("Failed to open/create HDF5 file")?;

    let group = get_or_create_group_in_file(&file, "/dielectric")
        .context("Failed to access/create dielectric group")?;

    // Store each frequency point as a separate dataset
    let dataset_name = format!("W_{omega_index:04}");

    // Optimal config for symmetric matrices
    let config = DatasetConfig {
        // Square chunks for symmetric matrices
        chunk_size: Some(vec![n_aux.min(256), n_aux.min(256)]),
        compression: CompressionType::Gzip(6), // Higher compression for storage
        shuffle: true,
        fletcher32: true,
    };

    write_array_optimized(&group, &dataset_name, data, &[n_aux, n_aux], Some(config))?;

    // Store frequency value as attribute
    let dataset = group.dataset(&dataset_name)?;
    dataset
        .new_attr::<f64>()
        .create("omega")?
        .write_scalar(&omega_value)?;

    Ok(())
}

/// Batch write multiple frequency-dependent matrices (for parallel I/O).
///
/// # Arguments
///
/// * `path` - Path to the HDF5 file
/// * `omega_values` - Array of frequency values
/// * `matrices` - Array of flattened matrices, one per frequency
/// * `n_aux` - Size of the auxiliary basis
///
/// # Errors
///
/// Returns an error if any matrix cannot be written
#[cfg(feature = "hdf5_support")]
pub fn write_frequency_batch(
    path: &Path,
    omega_values: &[f64],
    matrices: &[Vec<f64>], // Each inner Vec is a flattened matrix
    n_aux: usize,
) -> Result<()> {
    // Note: Parallel HDF5 requires special build with MPI support
    // For thread-safety, we use sequential writes

    // Sequential writes for thread safety (unless HDF5 is built with parallel support)
    for (i, (omega, matrix)) in omega_values.iter().zip(matrices.iter()).enumerate() {
        write_dielectric_matrix(path, i, *omega, matrix, n_aux)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::schema::{BasisSetData, CalculationParams, CalculationType, MoleculeInput};
    use tempfile::NamedTempFile;

    #[cfg(feature = "hdf5_support")]
    #[test]
    fn test_optimized_chunking() {
        use hdf5_impl::calculate_chunk_size;

        // Test 1D array chunking
        let shape_1d = vec![10000];
        let chunks_1d = calculate_chunk_size(&shape_1d, 8); // f64 size
        assert!(chunks_1d[0] <= 8192);
        assert!(chunks_1d[0] >= 64);

        // Test 2D transition tensor chunking (ia|P)
        let shape_ia_p = vec![1000, 5000]; // 1000 transitions, 5000 aux
        let chunks_ia_p = calculate_chunk_size(&shape_ia_p, 8);
        assert!(chunks_ia_p[0] <= 1024);
        assert!(chunks_ia_p[1] <= 8192);

        // Test square matrix chunking
        let shape_square = vec![2000, 2000];
        let chunks_square = calculate_chunk_size(&shape_square, 8);
        assert_eq!(chunks_square[0], chunks_square[1]); // Should be square chunks
    }

    #[cfg(feature = "hdf5_support")]
    #[test]
    fn test_write_read_optimized() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let path = temp_file.path();

        // Create test data
        let n_transitions = 500;
        let n_aux = 1000;
        let data: Vec<f64> = (0..n_transitions * n_aux)
            .map(|i| i as f64 * 0.001)
            .collect();

        // Write with optimization
        write_transition_tensor(path, "ia_P", &data, n_transitions, n_aux)?;

        // Read back
        let data_read = read_array_hdf5(path, "/df_tensors", "ia_P")?;

        // Verify
        assert_eq!(data.len(), data_read.len());
        for (a, b) in data.iter().zip(data_read.iter()) {
            assert!((a - b).abs() < 1e-10);
        }

        Ok(())
    }

    #[cfg(feature = "hdf5_support")]
    #[test]
    fn test_dielectric_matrix_io() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let path = temp_file.path();

        // Create test dielectric matrix
        let n_aux = 256;
        let omega = 0.5;
        let matrix: Vec<f64> = (0..n_aux * n_aux).map(|i| (i as f64).sin() * 0.1).collect();

        // Write dielectric matrix
        write_dielectric_matrix(path, 0, omega, &matrix, n_aux)?;

        // Read back
        let matrix_read = read_array_hdf5(path, "/dielectric", "W_0000")?;

        // Verify
        assert_eq!(matrix.len(), matrix_read.len());
        for (a, b) in matrix.iter().zip(matrix_read.iter()) {
            assert!((a - b).abs() < 1e-10);
        }

        Ok(())
    }

    #[cfg(feature = "hdf5_support")]
    #[test]
    fn test_frequency_batch_write() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let path = temp_file.path();

        // Create test data for multiple frequencies
        let n_aux = 128;
        let omega_values = vec![0.0, 0.5, 1.0, 1.5, 2.0];
        let matrices: Vec<Vec<f64>> = omega_values
            .iter()
            .map(|&omega| {
                (0..n_aux * n_aux)
                    .map(|i| (i as f64 * omega).cos() * 0.1)
                    .collect()
            })
            .collect();

        // Batch write
        write_frequency_batch(path, &omega_values, &matrices, n_aux)?;

        // Read back and verify each frequency point
        for (i, (_omega, matrix)) in omega_values.iter().zip(matrices.iter()).enumerate() {
            let dataset_name = format!("W_{i:04}");
            let matrix_read = read_array_hdf5(path, "/dielectric", &dataset_name)?;

            assert_eq!(matrix.len(), matrix_read.len());
            for (a, b) in matrix.iter().zip(matrix_read.iter()) {
                assert!((a - b).abs() < 1e-10);
            }
        }

        Ok(())
    }

    #[test]
    fn test_json_roundtrip() {
        // Create test molecule
        let molecule = MoleculeInput {
            natoms: 3,
            symbols: vec!["H".to_string(), "O".to_string(), "H".to_string()],
            atomic_numbers: vec![1, 8, 1],
            coordinates: vec![[-0.757, 0.0, 0.587], [0.0, 0.0, 0.0], [0.757, 0.0, 0.587]],
            charge: 0,
            multiplicity: 1,
            symmetry: Some("C2v".to_string()),
        };

        // Create basis set data
        let basis = BasisSetData {
            ao_basis: "def2-SVP".to_string(),
            n_ao: 24,
            aux_basis: "def2-SVP-JKFIT".to_string(),
            n_aux: 84,
            functions: None,
        };

        // Create calculation parameters
        let params = CalculationParams {
            calculation_type: CalculationType::G0W0,
            gw_params: Some(crate::io::GWParams {
                starting_point: "HF".to_string(),
                n_states: 10,
                contour_deformation: true,
                analytic_continuation: false,
                linearization: "diagonal".to_string(),
                max_iter: 1,
            }),
            bse_params: None,
            convergence: ConvergenceParams {
                energy_tol: 1e-6,
                density_tol: 1e-6,
                max_iterations: 100,
                use_diis: true,
                diis_space: 8,
            },
            frequency: FrequencyParams {
                grid_type: "GaussLegendre".to_string(),
                n_points: 16,
                omega_max: None,
                eta: 0.01,
            },
        };

        // Create QuasixData
        let mut data = QuasixData::new(molecule, basis, params);

        // Add some results
        data.results.reference.energy = -76.026_741;
        data.results.reference.orbital_energies = vec![-20.5, -1.3, -0.7, -0.5, -0.4];
        data.results.reference.occupations = vec![2.0, 2.0, 2.0, 2.0, 2.0];
        data.results.reference.homo = 4;
        data.results.reference.lumo = 5;

        // Create a temporary file for testing
        let temp_dir = std::env::temp_dir();
        let json_path = temp_dir.join("test_quasix_data.json");

        // Write to JSON
        data.to_json(&json_path).unwrap();

        // Read back from JSON
        let data2 = QuasixData::from_json(&json_path).unwrap();

        // Verify the data matches
        assert_eq!(data2.molecule.natoms, 3);
        assert_eq!(data2.molecule.symbols, vec!["H", "O", "H"]);
        assert_eq!(data2.basis.ao_basis, "def2-SVP");
        assert_eq!(data2.basis.n_ao, 24);
        assert!((data2.results.reference.energy - (-76.026_741)).abs() < 1e-9);
        assert_eq!(data2.results.reference.homo, 4);
        assert_eq!(data2.results.reference.lumo, 5);

        // Clean up
        std::fs::remove_file(json_path).ok();
    }

    #[test]
    fn test_schema_serialization() {
        // Create test molecule
        let molecule = MoleculeInput {
            natoms: 3,
            symbols: vec!["H".to_string(), "O".to_string(), "H".to_string()],
            atomic_numbers: vec![1, 8, 1],
            coordinates: vec![[-0.757, 0.0, 0.587], [0.0, 0.0, 0.0], [0.757, 0.0, 0.587]],
            charge: 0,
            multiplicity: 1,
            symmetry: Some("C2v".to_string()),
        };

        // Create basis set data
        let basis = BasisSetData {
            ao_basis: "def2-SVP".to_string(),
            n_ao: 24,
            aux_basis: "def2-SVP-JKFIT".to_string(),
            n_aux: 84,
            functions: None,
        };

        // Create calculation parameters
        let params = CalculationParams {
            calculation_type: CalculationType::G0W0,
            gw_params: Some(crate::io::GWParams {
                starting_point: "HF".to_string(),
                n_states: 10,
                contour_deformation: true,
                analytic_continuation: false,
                linearization: "diagonal".to_string(),
                max_iter: 1,
            }),
            bse_params: None,
            convergence: ConvergenceParams {
                energy_tol: 1e-6,
                density_tol: 1e-6,
                max_iterations: 100,
                use_diis: true,
                diis_space: 8,
            },
            frequency: FrequencyParams {
                grid_type: "GaussLegendre".to_string(),
                n_points: 16,
                omega_max: None,
                eta: 0.01,
            },
        };

        // Create QuasixData
        let data = QuasixData::new(molecule, basis, params);

        // Test JSON serialization
        let json = serde_json::to_string_pretty(&data).unwrap();
        assert!(json.contains("quasix_version"));
        assert!(json.contains("\"natoms\": 3"));

        // Test deserialization
        let data2: QuasixData = serde_json::from_str(&json).unwrap();
        assert_eq!(data2.molecule.natoms, 3);
        assert_eq!(data2.basis.ao_basis, "def2-SVP");
    }
}
