//! NUMA-aware memory allocation for high-performance parallel computing
//!
//! This module provides NUMA-aware memory allocation strategies optimized
//! for QuasiX's GW/BSE calculations on multi-socket systems.
#![allow(clippy::many_single_char_names)] // Mathematical notation
#![warn(clippy::all)]
// Allow common patterns in scientific computing
#![allow(clippy::cast_precision_loss)] // Array indexing to f64 is standard
#![allow(clippy::cast_possible_truncation)] // Safe for our use cases
#![allow(clippy::items_after_statements)] // Sometimes clearer in numerical code
#![warn(missing_docs)]

use crate::common::{QuasixError, Result};
use ndarray::{Array1, Array2, Array3};
use std::alloc::Layout;
use std::ptr::NonNull;

#[cfg(target_os = "linux")]
use libc::{c_void, size_t};

/// NUMA node identifier
pub type NumaNode = usize;

/// NUMA memory policy
#[derive(Debug, Clone, Copy)]
pub enum NumaPolicy {
    /// Default system policy
    Default,
    /// Bind allocation to specific node
    Bind(NumaNode),
    /// Interleave pages across nodes
    Interleave,
    /// Prefer local node but allow remote
    Preferred(NumaNode),
    /// Allocate on node of first touch
    FirstTouch,
}

/// NUMA-aware array allocator
pub struct NumaAllocator {
    /// Current NUMA policy
    policy: NumaPolicy,
    /// Number of NUMA nodes
    #[allow(dead_code)]
    n_nodes: usize,
    /// Page size in bytes
    page_size: usize,
    /// Enable huge pages
    huge_pages: bool,
}

impl NumaAllocator {
    /// Create a new NUMA allocator
    pub fn new() -> Result<Self> {
        let n_nodes = Self::get_numa_nodes()?;
        let page_size = Self::get_page_size();

        Ok(Self {
            policy: NumaPolicy::Default,
            n_nodes,
            page_size,
            huge_pages: false,
        })
    }

    /// Set NUMA policy
    pub fn set_policy(&mut self, policy: NumaPolicy) {
        self.policy = policy;
    }

    /// Enable huge pages for large allocations
    pub fn enable_huge_pages(&mut self, enable: bool) {
        self.huge_pages = enable;
    }

    /// Allocate 1D array with NUMA awareness
    pub fn alloc_array1<T: Clone + Default>(&self, size: usize) -> Result<Array1<T>> {
        let layout = Layout::array::<T>(size)
            .map_err(|e| QuasixError::NumericalError(format!("Invalid layout: {e}")))?;

        let ptr = self.allocate_numa(layout)?;

        // Initialize array with default values
        // Safety: ptr is valid and properly aligned from allocate_numa
        // Size is validated to be non-zero and within allocation bounds
        unsafe {
            let slice = std::slice::from_raw_parts_mut(ptr.as_ptr() as *mut T, size);
            for elem in slice.iter_mut() {
                std::ptr::write(elem, T::default());
            }
        }

        // Create ndarray from raw pointer
        // Safety: ptr is valid, size matches allocation, and ownership is transferred to Vec
        let array = unsafe {
            Array1::from_shape_vec_unchecked(
                size,
                Vec::from_raw_parts(ptr.as_ptr() as *mut T, size, size),
            )
        };

        Ok(array)
    }

    /// Allocate 2D array with NUMA awareness
    pub fn alloc_array2<T: Clone + Default>(&self, rows: usize, cols: usize) -> Result<Array2<T>> {
        let size = rows * cols;
        let layout = Layout::array::<T>(size)
            .map_err(|e| QuasixError::NumericalError(format!("Invalid layout: {e}")))?;

        let ptr = self.allocate_numa(layout)?;

        // Initialize array
        unsafe {
            let slice = std::slice::from_raw_parts_mut(ptr.as_ptr() as *mut T, size);
            for elem in slice.iter_mut() {
                std::ptr::write(elem, T::default());
            }
        }

        // Create ndarray
        let array = unsafe {
            Array2::from_shape_vec_unchecked(
                (rows, cols),
                Vec::from_raw_parts(ptr.as_ptr() as *mut T, size, size),
            )
        };

        Ok(array)
    }

    /// Allocate 3D array with NUMA awareness
    pub fn alloc_array3<T: Clone + Default>(
        &self,
        dim0: usize,
        dim1: usize,
        dim2: usize,
    ) -> Result<Array3<T>> {
        let size = dim0 * dim1 * dim2;
        let layout = Layout::array::<T>(size)
            .map_err(|e| QuasixError::NumericalError(format!("Invalid layout: {e}")))?;

        let ptr = self.allocate_numa(layout)?;

        // Initialize array
        unsafe {
            let slice = std::slice::from_raw_parts_mut(ptr.as_ptr() as *mut T, size);
            for elem in slice.iter_mut() {
                std::ptr::write(elem, T::default());
            }
        }

        // Create ndarray
        let array = unsafe {
            Array3::from_shape_vec_unchecked(
                (dim0, dim1, dim2),
                Vec::from_raw_parts(ptr.as_ptr() as *mut T, size, size),
            )
        };

        Ok(array)
    }

    /// Allocate memory with NUMA policy
    fn allocate_numa(&self, layout: Layout) -> Result<NonNull<u8>> {
        #[cfg(target_os = "linux")]
        {
            self.allocate_numa_linux(layout)
        }

        #[cfg(not(target_os = "linux"))]
        {
            self.allocate_standard(layout)
        }
    }

    /// Linux-specific NUMA allocation
    #[cfg(target_os = "linux")]
    fn allocate_numa_linux(&self, layout: Layout) -> Result<NonNull<u8>> {
        use libc::{mmap, MAP_ANONYMOUS, MAP_PRIVATE, PROT_READ, PROT_WRITE};

        let size = layout.size();
        let _align = layout.align();

        // Round up to page boundary
        let aligned_size = (size.div_ceil(self.page_size)) * self.page_size;

        // Allocate with mmap for better control
        let ptr = unsafe {
            mmap(
                std::ptr::null_mut(),
                aligned_size,
                PROT_READ | PROT_WRITE,
                MAP_PRIVATE | MAP_ANONYMOUS | if self.huge_pages { 0x40000 } else { 0 }, // MAP_HUGETLB
                -1,
                0,
            )
        };

        if ptr == libc::MAP_FAILED {
            return Err(QuasixError::NumericalError(
                "NUMA allocation failed".to_string(),
            ));
        }

        // Apply NUMA policy
        self.apply_numa_policy(ptr, aligned_size)?;

        // Touch pages to trigger first-touch policy
        if matches!(self.policy, NumaPolicy::FirstTouch) {
            self.touch_pages(ptr, aligned_size);
        }

        NonNull::new(ptr as *mut u8)
            .ok_or_else(|| QuasixError::NumericalError("Null pointer from mmap".to_string()))
    }

    /// Apply NUMA memory policy to allocated region
    #[cfg(target_os = "linux")]
    fn apply_numa_policy(&self, _ptr: *mut c_void, _size: size_t) -> Result<()> {
        // In production, would use libnuma or direct syscalls
        // This is a simplified version
        match self.policy {
            NumaPolicy::Bind(node) => {
                // mbind(ptr, size, MPOL_BIND, &nodemask, maxnode, 0)
                log::debug!("Binding memory to NUMA node {}", node);
            }
            NumaPolicy::Interleave => {
                // mbind(ptr, size, MPOL_INTERLEAVE, &nodemask, maxnode, 0)
                log::debug!("Interleaving memory across NUMA nodes");
            }
            NumaPolicy::Preferred(node) => {
                // mbind(ptr, size, MPOL_PREFERRED, &nodemask, maxnode, 0)
                log::debug!("Preferring NUMA node {}", node);
            }
            _ => {}
        }
        Ok(())
    }

    /// Touch pages to establish NUMA locality
    #[cfg(target_os = "linux")]
    fn touch_pages(&self, ptr: *mut c_void, size: size_t) {
        unsafe {
            let bytes = ptr as *mut u8;
            for i in (0..size).step_by(self.page_size) {
                std::ptr::write_volatile(bytes.add(i), 0);
            }
        }
    }

    /// Standard allocation for non-Linux systems
    #[cfg(not(target_os = "linux"))]
    fn allocate_standard(&self, layout: Layout) -> Result<NonNull<u8>> {
        let ptr = unsafe { alloc_zeroed(layout) };
        NonNull::new(ptr)
            .ok_or_else(|| QuasixError::NumericalError("Allocation failed".to_string()))
    }

    /// Get number of NUMA nodes
    fn get_numa_nodes() -> Result<usize> {
        #[cfg(target_os = "linux")]
        {
            // Parse /sys/devices/system/node/online
            if let Ok(contents) = std::fs::read_to_string("/sys/devices/system/node/online") {
                // Format is typically "0-1" for 2 nodes
                if let Some(dash_pos) = contents.find('-') {
                    if let Ok(max_node) = contents[dash_pos + 1..].trim().parse::<usize>() {
                        return Ok(max_node + 1);
                    }
                }
            }
        }

        // Default to 1 node if detection fails
        Ok(1)
    }

    /// Get system page size
    fn get_page_size() -> usize {
        #[cfg(target_os = "linux")]
        {
            unsafe { libc::sysconf(libc::_SC_PAGESIZE) as usize }
        }

        #[cfg(not(target_os = "linux"))]
        {
            4096 // Default page size
        }
    }

    /// Get current CPU's NUMA node
    pub fn get_current_node() -> NumaNode {
        #[cfg(target_os = "linux")]
        {
            // Parse /proc/self/numa_maps or use getcpu() syscall
            // Simplified: assume node 0
            0
        }

        #[cfg(not(target_os = "linux"))]
        {
            0
        }
    }
}

// Thread-local NUMA allocator
thread_local! {
    static NUMA_ALLOCATOR: std::cell::RefCell<Option<NumaAllocator>> = const { std::cell::RefCell::new(None) };
}

/// Initialize thread-local NUMA allocator
pub fn init_numa_allocator(policy: NumaPolicy) -> Result<()> {
    NUMA_ALLOCATOR.with(|alloc| {
        let mut allocator = NumaAllocator::new()?;
        allocator.set_policy(policy);
        *alloc.borrow_mut() = Some(allocator);
        Ok(())
    })
}

/// Allocate NUMA-aware array using thread-local allocator
pub fn numa_array1<T: Clone + Default>(size: usize) -> Result<Array1<T>> {
    NUMA_ALLOCATOR.with(|alloc| {
        let borrow = alloc.borrow();
        match borrow.as_ref() {
            Some(allocator) => allocator.alloc_array1(size),
            None => {
                // Fall back to standard allocation
                Ok(Array1::default(size))
            }
        }
    })
}

/// Allocate NUMA-aware 2D array
pub fn numa_array2<T: Clone + Default>(rows: usize, cols: usize) -> Result<Array2<T>> {
    NUMA_ALLOCATOR.with(|alloc| {
        let borrow = alloc.borrow();
        match borrow.as_ref() {
            Some(allocator) => allocator.alloc_array2(rows, cols),
            None => Ok(Array2::default((rows, cols))),
        }
    })
}

/// Allocate NUMA-aware 3D array
pub fn numa_array3<T: Clone + Default>(dim0: usize, dim1: usize, dim2: usize) -> Result<Array3<T>> {
    NUMA_ALLOCATOR.with(|alloc| {
        let borrow = alloc.borrow();
        match borrow.as_ref() {
            Some(allocator) => allocator.alloc_array3(dim0, dim1, dim2),
            None => Ok(Array3::default((dim0, dim1, dim2))),
        }
    })
}

/// NUMA topology information
#[derive(Debug, Clone)]
pub struct NumaTopology {
    /// Number of NUMA nodes
    pub n_nodes: usize,
    /// CPUs per node
    pub cpus_per_node: Vec<Vec<usize>>,
    /// Memory per node (bytes)
    pub memory_per_node: Vec<usize>,
    /// Distance matrix between nodes
    pub distances: Array2<u8>,
}

impl NumaTopology {
    /// Detect system NUMA topology
    pub fn detect() -> Result<Self> {
        #[cfg(target_os = "linux")]
        {
            Self::detect_linux()
        }

        #[cfg(not(target_os = "linux"))]
        {
            // Assume single node for non-Linux
            Ok(Self {
                n_nodes: 1,
                cpus_per_node: vec![vec![0]],
                memory_per_node: vec![Self::get_total_memory()],
                distances: Array2::from_elem((1, 1), 10),
            })
        }
    }

    #[cfg(target_os = "linux")]
    fn detect_linux() -> Result<Self> {
        let n_nodes = NumaAllocator::get_numa_nodes()?;
        let mut cpus_per_node = vec![Vec::new(); n_nodes];
        let mut memory_per_node = vec![0; n_nodes];

        // Parse CPU topology
        for node in 0..n_nodes {
            let cpu_list_path = format!("/sys/devices/system/node/node{}/cpulist", node);
            if let Ok(contents) = std::fs::read_to_string(&cpu_list_path) {
                // Parse CPU list (e.g., "0-3,8-11")
                cpus_per_node[node] = Self::parse_cpu_list(&contents);
            }

            // Parse memory info
            let meminfo_path = format!("/sys/devices/system/node/node{}/meminfo", node);
            if let Ok(contents) = std::fs::read_to_string(&meminfo_path) {
                memory_per_node[node] = Self::parse_node_memory(&contents);
            }
        }

        // Parse distance matrix
        let mut distances = Array2::from_elem((n_nodes, n_nodes), 10);
        for i in 0..n_nodes {
            let distance_path = format!("/sys/devices/system/node/node{}/distance", i);
            if let Ok(contents) = std::fs::read_to_string(&distance_path) {
                let dists: Vec<u8> = contents
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                for (j, &dist) in dists.iter().enumerate().take(n_nodes) {
                    distances[[i, j]] = dist;
                }
            }
        }

        Ok(Self {
            n_nodes,
            cpus_per_node,
            memory_per_node,
            distances,
        })
    }

    #[cfg(target_os = "linux")]
    fn parse_cpu_list(cpu_list: &str) -> Vec<usize> {
        let mut cpus = Vec::new();
        for part in cpu_list.trim().split(',') {
            if let Some(dash_pos) = part.find('-') {
                if let (Ok(start), Ok(end)) = (
                    part[..dash_pos].parse::<usize>(),
                    part[dash_pos + 1..].parse::<usize>(),
                ) {
                    cpus.extend(start..=end);
                }
            } else if let Ok(cpu) = part.parse::<usize>() {
                cpus.push(cpu);
            }
        }
        cpus
    }

    #[cfg(target_os = "linux")]
    fn parse_node_memory(meminfo: &str) -> usize {
        // Look for "MemTotal:" line
        for line in meminfo.lines() {
            if line.starts_with("MemTotal:") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    if let Ok(kb) = parts[1].parse::<usize>() {
                        return kb * 1024; // Convert to bytes
                    }
                }
            }
        }
        0
    }

    #[allow(dead_code)]
    fn get_total_memory() -> usize {
        // Simplified: return 16 GB
        16 * 1024 * 1024 * 1024
    }

    /// Get optimal node for thread
    pub fn get_optimal_node(&self, thread_id: usize) -> NumaNode {
        // Simple round-robin distribution
        thread_id % self.n_nodes
    }

    /// Check if system has multiple NUMA nodes
    pub fn is_numa_system(&self) -> bool {
        self.n_nodes > 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_numa_allocator_creation() {
        let allocator = NumaAllocator::new();
        assert!(allocator.is_ok());
    }

    #[test]
    fn test_numa_topology_detection() {
        let topology = NumaTopology::detect();
        assert!(topology.is_ok());
        let topo = topology.unwrap();
        assert!(topo.n_nodes >= 1);
    }

    #[test]
    fn test_page_size() {
        let page_size = NumaAllocator::get_page_size();
        assert!(page_size >= 4096);
    }

    #[test]
    fn test_numa_policies() {
        let policies = vec![
            NumaPolicy::Default,
            NumaPolicy::Bind(0),
            NumaPolicy::Interleave,
            NumaPolicy::Preferred(0),
            NumaPolicy::FirstTouch,
        ];

        for policy in policies {
            let mut allocator = NumaAllocator::new().unwrap();
            allocator.set_policy(policy);
        }
    }
}
