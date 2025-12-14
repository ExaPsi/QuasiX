"""
evGW convergence monitoring module for QuasiX.

This module provides real-time visualization and logging of evGW calculations,
with seamless PySCF integration and publication-quality plot generation.
"""

import numpy as np
from typing import Optional, Dict, Any, List, Tuple, Callable, Union
from dataclasses import dataclass, field
from pathlib import Path
from contextlib import contextmanager
import json
import time
import threading
import queue
import logging
import warnings
import tempfile
import shutil
import os
from datetime import datetime

# Configure logging
logger = logging.getLogger(__name__)

# Try importing matplotlib with graceful fallback
try:
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    import matplotlib.gridspec as gridspec
    from matplotlib.ticker import ScalarFormatter, LogFormatter
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    warnings.warn("Matplotlib not available. Plotting disabled, only logging active.")


@dataclass
class MonitorConfig:
    """Configuration for evGW monitoring."""
    # Plotting configuration
    plot_enabled: bool = True
    plot_update_interval: float = 0.5  # seconds
    plot_style: str = "publication"  # "publication", "interactive", "minimal"
    figure_size: Tuple[float, float] = (12, 8)
    dpi: int = 100
    save_plots: bool = False
    plot_dir: Path = field(default_factory=lambda: Path("./plots"))
    plot_format: str = "pdf"  # pdf, png, svg

    # Logging configuration
    log_enabled: bool = True
    log_file: Optional[Path] = None
    log_buffer_size: int = 10

    # Display configuration
    verbose: int = 1
    backend: str = "auto"  # matplotlib backend
    non_blocking: bool = True

    # Performance configuration
    memory_limit_mb: int = 100  # for log buffering
    use_memory_mapping: bool = False
    lazy_plotting: bool = True  # Batch plot updates

    # Convergence tracking
    track_oscillations: bool = True
    oscillation_window: int = 5

    def __post_init__(self):
        """Ensure paths are Path objects."""
        if self.plot_dir is not None:
            self.plot_dir = Path(self.plot_dir)
        if self.log_file is not None:
            self.log_file = Path(self.log_file)


@dataclass
class IterationUpdate:
    """Data from a single evGW iteration."""
    iteration: int
    timestamp: float
    qp_energies: np.ndarray
    z_factors: np.ndarray
    convergence: float  # Max energy change
    rms_convergence: float
    homo_energy: float
    lumo_energy: float
    gap: float
    timing: Dict[str, float]
    memory_mb: float
    damping: float = 0.5


class PlotManager:
    """Manages matplotlib backend and figure lifecycle."""

    def __init__(self, config: MonitorConfig):
        self.config = config
        self._original_backend = None
        self._figure = None
        self._axes = {}
        self._animation = None

    def setup_backend(self) -> str:
        """Configure matplotlib backend based on environment."""
        if not MATPLOTLIB_AVAILABLE:
            return 'none'

        if self.config.backend == "auto":
            backend = self._auto_select_backend()
        else:
            backend = self.config.backend

        self._original_backend = matplotlib.get_backend()
        try:
            matplotlib.use(backend)
        except Exception as e:
            logger.warning(f"Failed to set backend {backend}: {e}")
            backend = 'Agg'  # Fallback to non-interactive
            matplotlib.use(backend)

        return backend

    def _auto_select_backend(self) -> str:
        """Auto-detect best backend for current environment."""
        # Check if in Jupyter/IPython
        try:
            import IPython
            ipython = IPython.get_ipython()
            if ipython is not None:
                if 'IPKernelApp' in ipython.config:
                    return 'notebook'
        except ImportError:
            pass

        # Check display availability
        if not self._check_display():
            return 'Agg'

        # Try Qt, then Tk, then fallback to Agg
        for backend in ['Qt5Agg', 'TkAgg', 'Agg']:
            try:
                matplotlib.use(backend)
                return backend
            except:
                continue

        return 'Agg'

    def _check_display(self) -> bool:
        """Check if display is available."""
        if os.environ.get('DISPLAY'):
            return True
        if os.environ.get('WAYLAND_DISPLAY'):
            return True
        return False

    def create_dashboard(self) -> Tuple[plt.Figure, Dict[str, plt.Axes]]:
        """Create multi-panel dashboard for monitoring."""
        if not MATPLOTLIB_AVAILABLE:
            return None, {}

        # Apply publication style if requested
        if self.config.plot_style == "publication":
            self._apply_publication_style()

        # Create figure with constrained layout
        fig = plt.figure(figsize=self.config.figure_size,
                        dpi=self.config.dpi,
                        constrained_layout=True)

        # Create grid specification
        gs = gridspec.GridSpec(2, 2, figure=fig,
                              height_ratios=[1, 1],
                              width_ratios=[1, 1],
                              hspace=0.3, wspace=0.3)

        # Create subplots
        axes = {
            'convergence': fig.add_subplot(gs[0, :]),  # Top row, both columns
            'z_factors': fig.add_subplot(gs[1, 0]),    # Bottom left
            'gap': fig.add_subplot(gs[1, 1])           # Bottom right
        }

        # If we have more space, add oscillation detection panel
        if self.config.figure_size[0] > 12:
            # Recreate with 3x2 grid
            fig.clear()
            gs = gridspec.GridSpec(3, 2, figure=fig,
                                  height_ratios=[1, 1, 1],
                                  width_ratios=[1, 1],
                                  hspace=0.3, wspace=0.3)
            axes = {
                'convergence': fig.add_subplot(gs[0, :]),
                'z_factors': fig.add_subplot(gs[1, 0]),
                'gap': fig.add_subplot(gs[1, 1]),
                'oscillation': fig.add_subplot(gs[2, :])
            }

        self._setup_axes_properties(axes)
        self._figure = fig
        self._axes = axes

        return fig, axes

    def _setup_axes_properties(self, axes: Dict[str, plt.Axes]):
        """Configure axes for publication quality."""
        # Convergence plot
        if 'convergence' in axes:
            ax = axes['convergence']
            ax.set_xlabel('Iteration')
            ax.set_ylabel('Max |ΔE| (eV)')
            ax.set_yscale('log')
            ax.grid(True, alpha=0.3, which='both')
            ax.set_title('evGW Convergence', fontsize=12, fontweight='bold')

        # Z-factors plot
        if 'z_factors' in axes:
            ax = axes['z_factors']
            ax.set_xlabel('Iteration')
            ax.set_ylabel('Z Factor')
            ax.set_ylim(-0.05, 1.05)
            ax.axhline(y=0, color='r', linestyle='--', alpha=0.3, label='Physical bounds')
            ax.axhline(y=1, color='r', linestyle='--', alpha=0.3)
            ax.grid(True, alpha=0.3)
            ax.set_title('Quasiparticle Weights', fontsize=10)

        # Gap evolution plot
        if 'gap' in axes:
            ax = axes['gap']
            ax.set_xlabel('Iteration')
            ax.set_ylabel('HOMO-LUMO Gap (eV)')
            ax.grid(True, alpha=0.3)
            ax.set_title('Gap Evolution', fontsize=10)

        # Oscillation detection plot
        if 'oscillation' in axes:
            ax = axes['oscillation']
            ax.set_xlabel('Iteration')
            ax.set_ylabel('Oscillation Metric')
            ax.grid(True, alpha=0.3)
            ax.set_title('Convergence Pattern Analysis', fontsize=10)

    def _apply_publication_style(self):
        """Apply publication-quality settings."""
        plt.rcParams.update({
            'font.family': 'sans-serif',
            'font.sans-serif': ['Arial', 'DejaVu Sans'],
            'font.size': 10,
            'axes.titlesize': 12,
            'axes.labelsize': 11,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 9,
            'figure.titlesize': 14,
            'lines.linewidth': 2,
            'lines.markersize': 6,
            'axes.linewidth': 1,
            'grid.linewidth': 0.5,
            'xtick.major.width': 1,
            'ytick.major.width': 1
        })

    def cleanup(self):
        """Clean up matplotlib resources."""
        if self._animation is not None:
            self._animation.event_source.stop()
            self._animation = None
        if self._figure is not None:
            plt.close(self._figure)
            self._figure = None
        if self._original_backend is not None:
            matplotlib.use(self._original_backend)


class DataBuffer:
    """Efficient data buffering with numpy arrays."""

    def __init__(self, max_iterations: int = 1000, n_orbitals: int = 100):
        self.max_iterations = max_iterations
        self.n_orbitals = n_orbitals
        self._buffers = {}
        self._current_size = 0
        self._init_buffers()

    def _init_buffers(self):
        """Pre-allocate numpy arrays for efficiency."""
        # Scalar data per iteration
        self._buffers['iterations'] = np.zeros(self.max_iterations, dtype=np.int32)
        self._buffers['timestamps'] = np.zeros(self.max_iterations, dtype=np.float64)
        self._buffers['convergence'] = np.zeros(self.max_iterations, dtype=np.float64)
        self._buffers['rms_convergence'] = np.zeros(self.max_iterations, dtype=np.float64)
        self._buffers['gap'] = np.zeros(self.max_iterations, dtype=np.float64)
        self._buffers['homo_energy'] = np.zeros(self.max_iterations, dtype=np.float64)
        self._buffers['lumo_energy'] = np.zeros(self.max_iterations, dtype=np.float64)
        self._buffers['max_z'] = np.zeros(self.max_iterations, dtype=np.float64)
        self._buffers['min_z'] = np.zeros(self.max_iterations, dtype=np.float64)
        self._buffers['damping'] = np.zeros(self.max_iterations, dtype=np.float64)

        # Memory usage tracking
        self._buffers['memory_mb'] = np.zeros(self.max_iterations, dtype=np.float64)

        # For oscillation detection
        self._buffers['oscillation_metric'] = np.zeros(self.max_iterations, dtype=np.float64)

        # Array data per iteration (pre-allocate for max orbitals)
        self._buffers['qp_energies'] = np.zeros((self.max_iterations, self.n_orbitals), dtype=np.float64)
        self._buffers['z_factors'] = np.zeros((self.max_iterations, self.n_orbitals), dtype=np.float64)

    def add_iteration(self, update: IterationUpdate):
        """Add iteration data efficiently."""
        idx = self._current_size

        if idx >= self.max_iterations:
            self._expand_buffers()

        # Store scalar values
        self._buffers['iterations'][idx] = update.iteration
        self._buffers['timestamps'][idx] = update.timestamp
        self._buffers['convergence'][idx] = update.convergence
        self._buffers['rms_convergence'][idx] = update.rms_convergence
        self._buffers['gap'][idx] = update.gap
        self._buffers['homo_energy'][idx] = update.homo_energy
        self._buffers['lumo_energy'][idx] = update.lumo_energy
        self._buffers['max_z'][idx] = np.max(update.z_factors)
        self._buffers['min_z'][idx] = np.min(update.z_factors)
        self._buffers['damping'][idx] = update.damping
        self._buffers['memory_mb'][idx] = update.memory_mb

        # Store array data (pad or truncate as needed)
        n_qp = len(update.qp_energies) if update.qp_energies is not None else 0
        if n_qp > 0:
            n_store = min(n_qp, self.n_orbitals)
            self._buffers['qp_energies'][idx, :n_store] = update.qp_energies[:n_store]

        n_z = len(update.z_factors) if update.z_factors is not None else 0
        if n_z > 0:
            n_store = min(n_z, self.n_orbitals)
            self._buffers['z_factors'][idx, :n_store] = update.z_factors[:n_store]

        # Calculate oscillation metric if we have enough history
        if idx >= 2:
            self._calculate_oscillation_metric(idx)

        self._current_size += 1

    def _calculate_oscillation_metric(self, idx: int):
        """Calculate oscillation detection metric."""
        if idx < 2:
            return

        # Look at last 3 points for oscillation
        window = min(3, idx + 1)
        gaps = self._buffers['gap'][max(0, idx-window+1):idx+1]

        if len(gaps) >= 3:
            # Check for alternating pattern
            diffs = np.diff(gaps)
            if len(diffs) >= 2:
                # Oscillation if signs alternate
                oscillating = np.sign(diffs[0]) != np.sign(diffs[1])
                self._buffers['oscillation_metric'][idx] = 1.0 if oscillating else 0.0

    def _expand_buffers(self):
        """Expand buffers when full."""
        new_max = self.max_iterations * 2
        for key, arr in self._buffers.items():
            if arr.ndim == 1:
                new_arr = np.zeros(new_max, dtype=arr.dtype)
                new_arr[:self.max_iterations] = arr
            else:  # 2D array
                new_shape = (new_max,) + arr.shape[1:]
                new_arr = np.zeros(new_shape, dtype=arr.dtype)
                new_arr[:self.max_iterations] = arr
            self._buffers[key] = new_arr
        self.max_iterations = new_max

    def get_view(self, key: str) -> np.ndarray:
        """Return view of active data (no copy)."""
        return self._buffers[key][:self._current_size]

    def get_latest(self, key: str, n: int = 1):
        """Get last n values."""
        if self._current_size == 0:
            return np.array([])
        start = max(0, self._current_size - n)
        return self._buffers[key][start:self._current_size]


class JSONLogger:
    """Structured JSON logging for evGW monitoring."""

    def __init__(self, log_file: Path, buffer_size: int = 10):
        self.log_file = log_file
        self.buffer_size = buffer_size
        self._buffer = []
        self._metadata_written = False
        self._lock = threading.Lock()

    def write_metadata(self, metadata: Dict[str, Any]):
        """Write calculation metadata."""
        with self._lock:
            self._write_metadata_unlocked(metadata)

    def _write_metadata_unlocked(self, metadata: Dict[str, Any]):
        """Internal metadata write without locking."""
        with open(self.log_file, 'w') as f:
            json.dump({
                'metadata': metadata,
                'iterations': []
            }, f, indent=2)
        self._metadata_written = True

    def log_iteration(self, update: IterationUpdate):
        """Log iteration data with buffering."""
        iter_log = {
            'iteration': update.iteration,
            'timestamp': datetime.fromtimestamp(update.timestamp).isoformat(),
            'convergence': float(update.convergence),
            'rms_convergence': float(update.rms_convergence),
            'gap': float(update.gap),
            'homo_energy': float(update.homo_energy),
            'lumo_energy': float(update.lumo_energy),
            'min_z_factor': float(np.min(update.z_factors)),
            'max_z_factor': float(np.max(update.z_factors)),
            'damping': float(update.damping),
            'timing': update.timing,
            'memory_mb': float(update.memory_mb)
        }

        with self._lock:
            self._buffer.append(iter_log)

            if len(self._buffer) >= self.buffer_size:
                self._flush_unlocked()  # Call unlocked version since we already have the lock

    def _flush_unlocked(self):
        """Internal flush without locking - must be called while holding the lock."""
        if not self._buffer:
            return

        if not self._metadata_written:
            # Initialize file with empty structure
            with open(self.log_file, 'w') as f:
                json.dump({
                    'metadata': {'timestamp': datetime.now().isoformat()},
                    'iterations': []
                }, f, indent=2)
            self._metadata_written = True

        # Read existing data
        with open(self.log_file, 'r') as f:
            data = json.load(f)

        # Append buffered iterations
        data['iterations'].extend(self._buffer)

        # Write back
        with open(self.log_file, 'w') as f:
            json.dump(data, f, indent=2)

        self._buffer.clear()

    def flush(self):
        """Flush buffer to file."""
        with self._lock:
            self._flush_unlocked()

    def close(self):
        """Ensure all data is written."""
        self.flush()


class RealTimePlotter:
    """Handles real-time plot updates."""

    def __init__(self, plot_manager: PlotManager, data_buffer: DataBuffer):
        self.plot_manager = plot_manager
        self.data_buffer = data_buffer
        self.fig, self.axes = plot_manager.create_dashboard()
        self.plot_objects = {}
        self._init_plot_objects()

    def _init_plot_objects(self):
        """Initialize empty plot objects for updates."""
        if not MATPLOTLIB_AVAILABLE or self.fig is None:
            return

        # Convergence line
        if 'convergence' in self.axes:
            self.plot_objects['conv_line'], = self.axes['convergence'].plot(
                [], [], 'b-', linewidth=2, label='Max |ΔE|'
            )
            self.plot_objects['conv_rms'], = self.axes['convergence'].plot(
                [], [], 'g--', linewidth=1.5, label='RMS', alpha=0.7
            )

        # Z-factors lines
        if 'z_factors' in self.axes:
            self.plot_objects['z_min'], = self.axes['z_factors'].plot(
                [], [], 'b-', linewidth=2, label='Min Z', marker='o', markersize=4
            )
            self.plot_objects['z_max'], = self.axes['z_factors'].plot(
                [], [], 'r-', linewidth=2, label='Max Z', marker='s', markersize=4
            )

        # Gap line
        if 'gap' in self.axes:
            self.plot_objects['gap_line'], = self.axes['gap'].plot(
                [], [], 'g-', linewidth=2, marker='d', markersize=4
            )

        # Oscillation indicator
        if 'oscillation' in self.axes:
            self.plot_objects['osc_line'], = self.axes['oscillation'].plot(
                [], [], 'r-', linewidth=1.5, alpha=0.7
            )

    def update_plots(self):
        """Update all plots with current data."""
        if not MATPLOTLIB_AVAILABLE or self.fig is None:
            return

        # Get data views
        iterations = self.data_buffer.get_view('iterations')
        if len(iterations) == 0:
            return

        # Update convergence plot
        if 'convergence' in self.axes:
            self._update_convergence()

        # Update Z-factors plot
        if 'z_factors' in self.axes:
            self._update_z_factors()

        # Update gap plot
        if 'gap' in self.axes:
            self._update_gap()

        # Update oscillation plot
        if 'oscillation' in self.axes:
            self._update_oscillation()

        # Refresh canvas
        self.fig.canvas.draw_idle()

    def _update_convergence(self):
        """Update convergence plot."""
        ax = self.axes['convergence']
        iterations = self.data_buffer.get_view('iterations')
        convergence = self.data_buffer.get_view('convergence')
        rms_conv = self.data_buffer.get_view('rms_convergence')

        # Update main convergence line
        self.plot_objects['conv_line'].set_data(iterations, convergence * 27.2114)  # Convert to eV
        self.plot_objects['conv_rms'].set_data(iterations, rms_conv * 27.2114)

        # Auto-scale
        ax.relim()
        ax.autoscale_view()

        # Update legend
        if len(iterations) == 1:
            ax.legend(loc='best')

    def _update_z_factors(self):
        """Update Z-factors plot."""
        ax = self.axes['z_factors']
        iterations = self.data_buffer.get_view('iterations')
        min_z = self.data_buffer.get_view('min_z')
        max_z = self.data_buffer.get_view('max_z')

        self.plot_objects['z_min'].set_data(iterations, min_z)
        self.plot_objects['z_max'].set_data(iterations, max_z)

        # Auto-scale x-axis only
        ax.set_xlim(0, max(iterations[-1] + 1, 2))

        # Update legend
        if len(iterations) == 1:
            ax.legend(loc='best')

    def _update_gap(self):
        """Update gap evolution plot."""
        ax = self.axes['gap']
        iterations = self.data_buffer.get_view('iterations')
        gap = self.data_buffer.get_view('gap')

        self.plot_objects['gap_line'].set_data(iterations, gap * 27.2114)  # Convert to eV

        # Auto-scale
        ax.relim()
        ax.autoscale_view()

    def _update_oscillation(self):
        """Update oscillation detection plot."""
        if 'oscillation' not in self.axes:
            return

        ax = self.axes['oscillation']
        iterations = self.data_buffer.get_view('iterations')
        osc_metric = self.data_buffer.get_view('oscillation_metric')

        self.plot_objects['osc_line'].set_data(iterations, osc_metric)

        # Fixed y-axis for binary metric
        ax.set_ylim(-0.1, 1.1)
        ax.set_xlim(0, max(iterations[-1] + 1, 2))

    def save_figure(self, filename: str = None, dpi: int = 300):
        """Save current figure to file."""
        if not MATPLOTLIB_AVAILABLE or self.fig is None:
            return

        if filename is None:
            filename = f"evgw_monitor_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf"

        self.fig.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Figure saved to {filename}")


class EvGWMonitor:
    """
    Monitor and visualize evGW convergence in real-time.

    Integrates seamlessly with PySCF workflows and provides:
    - Real-time convergence plots
    - Z-factor evolution tracking
    - Energy level diagrams
    - JSON structured logging
    - Async/sync operation modes
    """

    def __init__(self, config: Optional[MonitorConfig] = None):
        self.config = config or MonitorConfig()
        self._is_running = False
        self._update_queue = queue.Queue()
        self._plot_thread = None
        self._callbacks = []
        self._finalized = False  # Track if finalize was called

        # Initialize components
        self.plot_manager = PlotManager(self.config) if self.config.plot_enabled else None
        self.data_buffer = DataBuffer()
        self.json_logger = JSONLogger(self.config.log_file) if self.config.log_file else None
        self.real_time_plotter = None

        # Setup backend if plotting enabled
        if self.config.plot_enabled and MATPLOTLIB_AVAILABLE:
            backend = self.plot_manager.setup_backend()
            logger.info(f"Matplotlib backend: {backend}")

    def __del__(self):
        """Destructor to ensure cleanup even if not properly finalized."""
        # Only cleanup if not already finalized
        if not getattr(self, '_finalized', True):  # Default to True if attribute missing
            try:
                # Emergency stop without saving
                if getattr(self, '_is_running', False):
                    self._is_running = False
                    # Try to stop thread
                    if hasattr(self, '_plot_thread') and self._plot_thread and self._plot_thread.is_alive():
                        if hasattr(self, '_update_queue'):
                            self._update_queue.put(None)
                        # Don't wait too long in destructor
                        self._plot_thread.join(timeout=0.5)
                    # Cleanup plot manager
                    if hasattr(self, 'plot_manager') and self.plot_manager:
                        self.plot_manager.cleanup()
            except Exception:
                pass  # Ignore errors in destructor

    def attach_to_evgw(self, evgw_instance):
        """Attach monitor to an EvGW calculation instance."""
        evgw_instance.monitor = self
        # Register our update method as callback
        if hasattr(evgw_instance, 'register_callback'):
            evgw_instance.register_callback(self.update)

    def register_callback(self, callback: Callable):
        """Register additional callback for iteration updates."""
        self._callbacks.append(callback)

    def update(self, iteration_data: Dict[str, Any]):
        """
        Update monitor with new iteration data.

        Parameters
        ----------
        iteration_data : dict
            Dictionary containing iteration information
        """
        # Convert to IterationUpdate object
        update = self._create_update_object(iteration_data)

        # Add to buffer
        self.data_buffer.add_iteration(update)

        # Log to JSON if enabled
        if self.json_logger:
            self.json_logger.log_iteration(update)

        # Queue for plotting if non-blocking
        if self.config.non_blocking and self.config.plot_enabled:
            self._update_queue.put(update)
        elif self.config.plot_enabled and self.real_time_plotter:
            self.real_time_plotter.update_plots()

        # Call registered callbacks
        for callback in self._callbacks:
            try:
                callback(update)
            except Exception as e:
                logger.warning(f"Callback failed: {e}")

    def _create_update_object(self, data: Dict[str, Any]) -> IterationUpdate:
        """Convert raw iteration data to IterationUpdate object."""
        # Extract required fields with defaults
        return IterationUpdate(
            iteration=data.get('iteration', data.get('cycle', 0)),
            timestamp=time.time(),
            qp_energies=np.asarray(data.get('qp_energies', [])),
            z_factors=np.asarray(data.get('z_factors', [])),
            convergence=float(data.get('convergence', data.get('energy_change', 0.0))),
            rms_convergence=float(data.get('rms_change', 0.0)),
            homo_energy=float(data.get('homo_energy', 0.0)),
            lumo_energy=float(data.get('lumo_energy', 0.0)),
            gap=float(data.get('gap', 0.0)),
            timing=data.get('timing', {}),
            memory_mb=float(data.get('memory_mb', 0.0)),
            damping=float(data.get('damping_used', 0.5))
        )

    @contextmanager
    def monitoring_context(self):
        """Context manager for automatic setup/cleanup."""
        try:
            self.start()
            yield self
        finally:
            self.stop()

    def start(self):
        """Start monitoring (initialize plots, start threads)."""
        if self._is_running:
            return

        self._is_running = True

        # Initialize plotter if enabled
        if self.config.plot_enabled and MATPLOTLIB_AVAILABLE:
            self.real_time_plotter = RealTimePlotter(self.plot_manager, self.data_buffer)

            # Start plot update thread if non-blocking
            if self.config.non_blocking:
                self._plot_thread = threading.Thread(target=self._plot_update_loop)
                self._plot_thread.daemon = True
                self._plot_thread.start()

        # Write metadata if logging
        if self.json_logger:
            metadata = {
                'timestamp': datetime.now().isoformat(),
                'config': {
                    'plot_enabled': self.config.plot_enabled,
                    'update_interval': self.config.plot_update_interval,
                    'backend': self.plot_manager.setup_backend() if self.plot_manager else 'none'
                }
            }
            self.json_logger.write_metadata(metadata)

        logger.info("evGW monitoring started")

    def stop(self):
        """Stop monitoring and clean up resources."""
        if not self._is_running:
            return

        self._is_running = False

        # Stop plot thread first with proper cleanup
        if self._plot_thread and self._plot_thread.is_alive():
            # Signal thread to stop by setting flag and sending stop signal
            self._update_queue.put(None)  # Send stop signal
            # Give thread time to finish gracefully
            self._plot_thread.join(timeout=2.0)
            # If thread is still alive, it's stuck - clear the queue
            if self._plot_thread.is_alive():
                # Clear any remaining items in queue
                while not self._update_queue.empty():
                    try:
                        self._update_queue.get_nowait()
                    except queue.Empty:
                        break
                logger.warning("Plot thread did not stop gracefully")

        # Clear the thread reference
        self._plot_thread = None

        logger.info("evGW monitoring stopped")

    def _plot_update_loop(self):
        """Background thread for plot updates."""
        last_update = time.time()

        while self._is_running:
            try:
                # Get update from queue with short timeout to check _is_running frequently
                update = self._update_queue.get(timeout=0.1)

                if update is None:  # Stop signal
                    break

                # Double-check we're still running after getting an update
                if not self._is_running:
                    break

                # Check if enough time has passed
                current_time = time.time()
                if current_time - last_update >= self.config.plot_update_interval:
                    if self.real_time_plotter and self._is_running:
                        try:
                            self.real_time_plotter.update_plots()
                            last_update = current_time
                        except Exception as e:
                            logger.warning(f"Plot update failed: {e}")
                            # Don't crash the thread on plot errors

            except queue.Empty:
                # Check if we should stop even without updates
                if not self._is_running:
                    break
                continue
            except Exception as e:
                logger.warning(f"Unexpected error in plot thread: {e}")
                if not self._is_running:
                    break

    def get_convergence_data(self) -> Dict[str, np.ndarray]:
        """Get convergence history data."""
        return {
            'iterations': self.data_buffer.get_view('iterations'),
            'convergence': self.data_buffer.get_view('convergence'),
            'rms_convergence': self.data_buffer.get_view('rms_convergence'),
            'gap': self.data_buffer.get_view('gap'),
            'z_min': self.data_buffer.get_view('min_z'),
            'z_max': self.data_buffer.get_view('max_z')
        }

    def save_data(self, filename: str):
        """Save monitoring data to file."""
        data = self.get_convergence_data()

        if filename.endswith('.npz'):
            np.savez(filename, **data)
        elif filename.endswith('.json'):
            # Convert to JSON-serializable format
            json_data = {k: v.tolist() for k, v in data.items()}
            with open(filename, 'w') as f:
                json.dump(json_data, f, indent=2)
        else:
            # Default to numpy format
            np.savez(filename + '.npz', **data)

        logger.info(f"Monitoring data saved to {filename}")

    # === Additional methods for S5-3 testing compatibility ===

    def update_iteration(self, update: IterationUpdate):
        """Add an iteration update (alias for update method)."""
        # Convert IterationUpdate to dict format expected by update()
        data = {
            'iteration': update.iteration,
            'qp_energies': update.qp_energies,
            'z_factors': update.z_factors,
            'convergence': update.convergence,
            'rms_change': update.rms_convergence,
            'homo_energy': update.homo_energy,
            'lumo_energy': update.lumo_energy,
            'gap': update.gap,
            'timing': update.timing,
            'memory_mb': update.memory_mb,
            'damping_used': update.damping
        }
        self.update(data)

    @property
    def n_iterations(self) -> int:
        """Get number of iterations completed."""
        return self.data_buffer._current_size

    @property
    def iterations(self) -> List[IterationUpdate]:
        """Get list of all iteration updates."""
        # Reconstruct IterationUpdate objects from buffer
        iterations = []
        for i in range(self.n_iterations):
            update = IterationUpdate(
                iteration=int(self.data_buffer._buffers['iterations'][i]),
                timestamp=self.data_buffer._buffers['timestamps'][i],
                qp_energies=self.data_buffer._buffers['qp_energies'][i],
                z_factors=self.data_buffer._buffers['z_factors'][i],
                convergence=self.data_buffer._buffers['convergence'][i],
                rms_convergence=self.data_buffer._buffers['rms_convergence'][i],
                homo_energy=self.data_buffer._buffers['homo_energy'][i],
                lumo_energy=self.data_buffer._buffers['lumo_energy'][i],
                gap=self.data_buffer._buffers['gap'][i],
                timing={},
                memory_mb=0.0
            )
            iterations.append(update)
        return iterations

    def get_convergence_metric(self, iteration: int) -> float:
        """Get convergence metric for specific iteration."""
        if iteration >= self.n_iterations:
            return 0.0
        return self.data_buffer._buffers['convergence'][iteration]

    def is_converged(self, threshold: float = 1e-4) -> bool:
        """Check if calculation has converged."""
        if self.n_iterations < 2:
            return False
        return self.data_buffer._buffers['convergence'][self.n_iterations - 1] < threshold

    def get_gap(self, iteration: int) -> float:
        """Get HOMO-LUMO gap at specific iteration."""
        if iteration >= self.n_iterations:
            return 0.0
        return self.data_buffer._buffers['gap'][iteration]

    def get_rms_convergence(self, iteration: int) -> float:
        """Get RMS convergence at specific iteration."""
        if iteration >= self.n_iterations:
            return 0.0
        return self.data_buffer._buffers['rms_convergence'][iteration]

    def check_z_factor_bounds(self, dsigma: float) -> List[str]:
        """Check Z-factor physical bounds."""
        warnings_list = []
        if dsigma >= 1.0:
            warnings_list.append("Z-factor would be negative (unphysical)")
        elif dsigma <= -1.0:
            warnings_list.append("Z-factor would exceed 1 (unphysical)")
        elif dsigma < -0.99:
            z = 1.0 / (1.0 - dsigma)
            warnings_list.append(f"Z-factor very small ({z:.3f}), near-metallic behavior")
        return warnings_list

    def compute_qp_energy(self, epsilon: float, sigma: float, vxc: float, z: float) -> float:
        """Compute quasiparticle energy using Z-factor."""
        return epsilon + z * (sigma - vxc)

    def detect_oscillation(self) -> Tuple[bool, Optional[int]]:
        """Detect oscillations in convergence."""
        if self.n_iterations < 5:
            return False, None

        # Get recent HOMO energies
        homo_data = self.data_buffer.get_view('homo_energy')

        # Simple oscillation detection: check for alternating pattern
        diffs = np.diff(homo_data)
        sign_changes = np.diff(np.sign(diffs))
        n_sign_changes = np.sum(sign_changes != 0)

        # If 70% or more are sign changes, likely oscillating
        if n_sign_changes >= 0.7 * len(sign_changes):
            # Try to detect period
            for period in [2, 3, 4, 5]:
                if self._check_period(homo_data, period):
                    return True, period
            return True, None

        return False, None

    def _check_period(self, data: np.ndarray, period: int) -> bool:
        """Check if data has specific period."""
        if len(data) < 2 * period:
            return False

        # Check if pattern repeats with given period
        # Compare the last 'period' elements with the ones before
        n = len(data)

        # Need at least 2 full periods to check
        if n < 2 * period:
            return False

        # Check if the pattern repeats for the last few cycles
        n_cycles = n // period
        if n_cycles < 2:
            return False

        # Extract last few complete cycles
        last_cycles = data[-(n_cycles * period):]

        # Reshape to compare cycles
        cycles = last_cycles.reshape(n_cycles, period)

        # Check if all cycles are similar (within tolerance)
        first_cycle = cycles[0]
        for i in range(1, n_cycles):
            if not np.allclose(cycles[i], first_cycle, rtol=0.15, atol=0.05):
                return False

        return True

    def compute_autocorrelation(self, signal: np.ndarray) -> np.ndarray:
        """Compute autocorrelation of signal."""
        n = len(signal)
        signal = signal - np.mean(signal)
        autocorr = np.correlate(signal, signal, mode='full')
        autocorr = autocorr[n-1:]
        autocorr = autocorr / autocorr[0]  # Normalize
        return autocorr

    def find_autocorrelation_peaks(self, autocorr: np.ndarray) -> List[int]:
        """Find peaks in autocorrelation."""
        from scipy.signal import find_peaks
        peaks, _ = find_peaks(autocorr[1:], height=0.5)  # Skip lag 0
        return [p + 1 for p in peaks]  # Adjust for skipped lag 0

    def compute_window_variance(self, data: List[float], start: int, window_size: int) -> float:
        """Compute variance over sliding window."""
        end = min(start + window_size, len(data))
        if end - start < 2:
            return 0.0
        window_data = data[start:end]
        return np.var(window_data)

    def estimate_convergence_rate(self) -> float:
        """Estimate exponential convergence rate."""
        if self.n_iterations < 3:
            return 0.0

        convergence = self.data_buffer.get_view('convergence')
        # Take log of convergence values (skip zeros)
        nonzero = convergence > 1e-10
        if np.sum(nonzero) < 2:
            return 0.0

        log_conv = np.log(convergence[nonzero])
        iterations = np.arange(len(convergence))[nonzero]

        # Fit linear regression to log(error) vs iteration
        from scipy.stats import linregress
        slope, _, _, _, _ = linregress(iterations, log_conv)

        # Rate = exp(slope)
        return np.exp(slope)

    def export_json(self) -> Dict[str, Any]:
        """Export monitoring data as JSON dict."""
        return {
            'metadata': {
                'version': '1.0',
                'timestamp': datetime.now().isoformat(),
                'system': 'QuasiX evGW'
            },
            'iterations': [
                {
                    'iteration': int(i),
                    'homo': float(self.data_buffer._buffers['homo_energy'][i]),
                    'lumo': float(self.data_buffer._buffers['lumo_energy'][i]),
                    'gap': float(self.data_buffer._buffers['gap'][i]),
                    'z_factors': self.data_buffer._buffers['z_factors'][i].tolist() if hasattr(self.data_buffer._buffers['z_factors'][i], 'tolist') else [],
                    'convergence': float(self.data_buffer._buffers['convergence'][i]),
                    'rms_convergence': float(self.data_buffer._buffers['rms_convergence'][i])
                }
                for i in range(self.n_iterations)
            ],
            'summary': {
                'converged': bool(self.is_converged()),  # Convert np.bool_ to native bool
                'n_iterations': int(self.n_iterations),
                'final_gap': float(self.data_buffer._buffers['gap'][self.n_iterations - 1]) if self.n_iterations > 0 else 0.0
            }
        }

    def export_json_string(self) -> str:
        """Export monitoring data as JSON string."""
        return json.dumps(self.export_json(), indent=2)

    @classmethod
    def from_json(cls, json_str: str) -> 'EvGWMonitor':
        """Create monitor from JSON string."""
        data = json.loads(json_str)
        monitor = cls()

        # Populate from JSON data
        for it_data in data['iterations']:
            update = IterationUpdate(
                iteration=it_data['iteration'],
                timestamp=0.0,
                qp_energies=np.zeros(2),  # Placeholder
                z_factors=np.array(it_data.get('z_factors', [])),
                convergence=it_data.get('convergence', 0.0),
                rms_convergence=it_data.get('rms_convergence', 0.0),
                homo_energy=it_data['homo'],
                lumo_energy=it_data['lumo'],
                gap=it_data['gap'],
                timing={},
                memory_mb=0.0
            )
            monitor.update_iteration(update)

        return monitor

    def get_homo_history(self) -> np.ndarray:
        """Get HOMO energy history."""
        if self.n_iterations == 0:
            return np.array([])
        return self.data_buffer.get_view('homo_energy')

    def get_lumo_history(self) -> np.ndarray:
        """Get LUMO energy history."""
        if self.n_iterations == 0:
            return np.array([])
        return self.data_buffer.get_view('lumo_energy')

    def mark_converged(self):
        """Mark calculation as converged."""
        self._converged = True

    def mark_failed(self, reason: str):
        """Mark calculation as failed."""
        self._converged = False
        self._failure_reason = reason

    def get_failure_reason(self) -> str:
        """Get failure reason if failed."""
        return getattr(self, '_failure_reason', '')

    def get_recommendations(self) -> str:
        """Get recommendations based on convergence behavior."""
        recommendations = []

        # Check for oscillations
        is_osc, period = self.detect_oscillation()
        if is_osc:
            recommendations.append(f"Oscillations detected (period={period}). Consider increasing damping.")

        # Check Z-factors
        if self.n_iterations > 0:
            min_z = self.data_buffer._buffers['min_z'][self.n_iterations - 1]
            if min_z < 0.5:
                recommendations.append(f"Small Z-factors ({min_z:.2f}) indicate near-metallic behavior.")

        # Check convergence rate
        if self.n_iterations > 10 and not self.is_converged():
            recommendations.append("Slow convergence. Consider adjusting damping or initial guess.")

        return " ".join(recommendations)

    def get_warnings(self) -> List[str]:
        """Get list of warnings."""
        return getattr(self, '_warnings', [])

    def generate_convergence_plot(self) -> Optional[plt.Figure]:
        """Generate convergence plot."""
        if not MATPLOTLIB_AVAILABLE:
            return None

        try:
            fig, axes = self.plot_manager.create_dashboard() if self.plot_manager else (None, {})
            if fig is None:
                return None

            # Update with current data
            if self.real_time_plotter:
                self.real_time_plotter.update_plots()

            return fig
        except Exception as e:
            logger.warning(f"Failed to generate plot: {e}")
            return None

    def estimate_memory_usage(self) -> float:
        """Estimate memory usage in MB."""
        # Rough estimate based on buffer sizes
        n_elements = self.n_iterations * (self.data_buffer.n_orbitals + 10)
        bytes_per_element = 8  # float64
        return (n_elements * bytes_per_element) / (1024 * 1024)

    def has_oscillations(self) -> bool:
        """Check if oscillations have been detected.

        Returns:
            True if oscillations detected
        """
        is_osc, _ = self.detect_oscillation()
        return is_osc

    def get_convergence_metrics(self) -> Dict[str, Any]:
        """Get comprehensive convergence metrics.

        Returns:
            Dictionary with convergence statistics
        """
        if self.n_iterations == 0:
            return {
                'converged': False,
                'n_iterations': 0,
                'oscillations_detected': False,
                'final_error': 0.0,
                'final_rms_error': 0.0,
                'runtime': 0.0
            }

        # Get final values
        final_idx = self.n_iterations - 1
        convergence = self.data_buffer._buffers['convergence'][final_idx]
        rms_convergence = self.data_buffer._buffers['rms_convergence'][final_idx]

        # Check oscillations
        is_osc, period = self.detect_oscillation()

        # Calculate runtime
        if hasattr(self, '_start_time'):
            runtime = time.time() - self._start_time
        else:
            runtime = 0.0

        # Get gap evolution
        gap_initial = self.data_buffer._buffers['gap'][0] if self.n_iterations > 0 else 0.0
        gap_final = self.data_buffer._buffers['gap'][final_idx]
        gap_change = gap_final - gap_initial

        # Get Z-factor statistics
        z_min = self.data_buffer._buffers['min_z'][final_idx] if 'min_z' in self.data_buffer._buffers else 0.0
        z_max = self.data_buffer._buffers['max_z'][final_idx] if 'max_z' in self.data_buffer._buffers else 1.0

        # Convergence rate
        conv_rate = self.estimate_convergence_rate()

        metrics = {
            'converged': bool(self.is_converged()),  # Convert np.bool_ to native bool
            'n_iterations': int(self.n_iterations),
            'oscillations_detected': bool(is_osc),  # Convert np.bool_ to native bool
            'oscillation_period': int(period) if period is not None else None,
            'final_error': float(convergence),
            'final_rms_error': float(rms_convergence),
            'runtime': float(runtime),
            'gap_initial': float(gap_initial),
            'gap_final': float(gap_final),
            'gap_change': float(gap_change),
            'z_min': float(z_min),
            'z_max': float(z_max),
            'convergence_rate': float(conv_rate),
            'memory_mb': float(self.estimate_memory_usage())
        }

        return metrics

    def finalize(self, save_plots: bool = True, save_json: bool = True) -> Dict[str, Any]:
        """Finalize monitoring and generate summary.

        Args:
            save_plots: Save plots to file
            save_json: Save JSON log

        Returns:
            Final summary dictionary
        """
        # Mark as finalized to prevent duplicate cleanup
        if self._finalized:
            logger.debug("Monitor already finalized, skipping")
            return self.get_convergence_metrics()

        self._finalized = True

        # Stop monitoring first to ensure threads are properly shut down
        self.stop()

        # Get final metrics
        summary = self.get_convergence_metrics()

        # Add recommendations
        summary['recommendations'] = self.get_recommendations()
        summary['warnings'] = self.get_warnings()

        # Save plots if requested (after stop() has been called)
        if save_plots and self.config.save_plots:
            if self.real_time_plotter:
                try:
                    plot_file = self.config.plot_dir / f"evgw_final_{datetime.now().strftime('%Y%m%d_%H%M%S')}.{self.config.plot_format}"
                    self.config.plot_dir.mkdir(parents=True, exist_ok=True)
                    self.real_time_plotter.save_figure(str(plot_file))
                    summary['plot_file'] = str(plot_file)
                    logger.info(f"Final plot saved to {plot_file}")
                except Exception as e:
                    logger.warning(f"Failed to save final plot in finalize: {e}")

        # Save JSON log if requested (after stop() has been called)
        if save_json and self.json_logger:
            try:
                # Add summary to JSON log
                self.json_logger.write_metadata({'summary': summary})
                self.json_logger.close()
                summary['json_file'] = str(self.config.log_file)
                logger.info(f"JSON log saved to {self.config.log_file}")
            except Exception as e:
                logger.warning(f"Failed to save JSON log in finalize: {e}")

        # Clean up matplotlib resources
        if self.plot_manager:
            try:
                self.plot_manager.cleanup()
            except Exception as e:
                logger.warning(f"Failed to cleanup plot manager: {e}")

        # Log summary if verbose
        if self.config.verbose > 0:
            logger.info("=" * 60)
            logger.info("evGW Monitoring Summary:")
            logger.info(f"  Iterations: {summary['n_iterations']}")
            logger.info(f"  Runtime: {summary['runtime']:.2f} s")
            logger.info(f"  Converged: {summary['converged']}")
            logger.info(f"  Final RMS: {summary['final_rms_error']:.2e} Ha")
            logger.info(f"  Gap change: {summary['gap_change']:.3f} Ha")
            logger.info(f"  Oscillations: {summary['oscillations_detected']}")
            if summary['recommendations']:
                logger.info(f"  Recommendations: {summary['recommendations']}")
            logger.info("=" * 60)

        return summary

    def log_iteration(self, iteration: int, qp_energies: np.ndarray,
                     z_factors: np.ndarray, qp_energies_prev: Optional[np.ndarray] = None,
                     damping: float = 1.0, converged: bool = False,
                     homo_idx: Optional[int] = None, lumo_idx: Optional[int] = None) -> Dict[str, Any]:
        """Log iteration data with minimal overhead.

        Args:
            iteration: Iteration number
            qp_energies: Current QP energies
            z_factors: Current Z factors
            qp_energies_prev: Previous QP energies (for computing changes)
            damping: Damping factor used
            converged: Whether iteration converged
            homo_idx: HOMO index (optional, will be guessed if not provided)
            lumo_idx: LUMO index (optional, will be guessed if not provided)

        Returns:
            Dictionary with convergence metrics
        """
        t_start = time.perf_counter()

        # Guess HOMO/LUMO if not provided
        if homo_idx is None:
            homo_idx = len(qp_energies) // 2 - 1
        if lumo_idx is None:
            lumo_idx = len(qp_energies) // 2

        # Compute changes
        if qp_energies_prev is not None:
            changes = qp_energies - qp_energies_prev
            max_change = np.max(np.abs(changes))
            rms_change = np.sqrt(np.mean(changes**2))
        else:
            max_change = rms_change = 0.0

        # Create iteration update
        update = IterationUpdate(
            iteration=iteration,
            timestamp=time.time(),
            qp_energies=qp_energies,
            z_factors=z_factors,
            convergence=max_change,
            rms_convergence=rms_change,
            homo_energy=qp_energies[homo_idx],
            lumo_energy=qp_energies[lumo_idx],
            gap=qp_energies[lumo_idx] - qp_energies[homo_idx],
            timing={'iteration': time.perf_counter() - t_start},
            memory_mb=self.estimate_memory_usage(),
            damping=damping
        )

        # Update through existing method
        self.update_iteration(update)

        # Performance check
        dt = (time.perf_counter() - t_start) * 1e6  # microseconds
        if dt > 100:
            logger.debug(f"Monitor update took {dt:.0f} μs (target: <100 μs)")

        return {
            'max_change': max_change,
            'rms_change': rms_change,
            'homo_lumo_gap': qp_energies[lumo_idx] - qp_energies[homo_idx],
            'mean_z': np.mean(z_factors),
            'oscillations': self.has_oscillations()
        }


class LogAnalyzer:
    """Utilities for analyzing evGW monitoring logs."""

    def __init__(self, log_file: Path):
        self.log_file = Path(log_file)
        self.data = self._load_log()

    def _load_log(self) -> Dict[str, Any]:
        """Load and validate log file."""
        with open(self.log_file, 'r') as f:
            data = json.load(f)
        return data

    def get_convergence_history(self) -> Tuple[np.ndarray, np.ndarray]:
        """Extract convergence history."""
        iterations = [it['iteration'] for it in self.data['iterations']]
        convergence = [it['convergence'] for it in self.data['iterations']]
        return np.array(iterations), np.array(convergence)

    def analyze_z_factors(self) -> Dict[str, Any]:
        """Analyze Z-factor evolution."""
        z_min = [it['min_z_factor'] for it in self.data['iterations']]
        z_max = [it['max_z_factor'] for it in self.data['iterations']]

        return {
            'min_trajectory': np.array(z_min),
            'max_trajectory': np.array(z_max),
            'final_range': (z_min[-1], z_max[-1]) if z_min else (0, 1),
            'converged': all(0 <= z <= 1 for z in z_min + z_max)
        }

    def export_for_publication(self, output_dir: Path):
        """Export data in publication-ready format."""
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)

        # Export convergence data
        iterations, convergence = self.get_convergence_history()
        np.savetxt(output_dir / 'convergence.dat',
                  np.column_stack([iterations, convergence]),
                  header='Iteration Convergence(Ha)',
                  fmt='%d %.6e')

        # Export Z-factor analysis
        z_analysis = self.analyze_z_factors()
        np.savetxt(output_dir / 'z_factors.dat',
                  np.column_stack([z_analysis['min_trajectory'],
                                 z_analysis['max_trajectory']]),
                  header='Min_Z Max_Z',
                  fmt='%.6f %.6f')

        logger.info(f"Publication data exported to {output_dir}")


def create_publication_plots(analyzer: LogAnalyzer, save_dir: Path = None) -> plt.Figure:
    """Create publication-quality plots from log data."""
    if not MATPLOTLIB_AVAILABLE:
        logger.warning("Matplotlib not available, cannot create plots")
        return None

    # Set publication style
    plt.style.use('seaborn-v0_8-paper')
    plt.rcParams.update({
        'font.size': 10,
        'axes.labelsize': 11,
        'axes.titlesize': 12,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'figure.titlesize': 14
    })

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)

    # Convergence plot
    iterations, convergence = analyzer.get_convergence_history()
    ax = axes[0, 0]
    ax.semilogy(iterations, convergence * 27.2114, 'b-', linewidth=2)
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Max |ΔE| (eV)')
    ax.set_title('evGW Convergence')
    ax.grid(True, alpha=0.3)

    # Z-factors plot
    z_analysis = analyzer.analyze_z_factors()
    ax = axes[0, 1]
    ax.plot(iterations, z_analysis['min_trajectory'], 'b-', label='Min Z', linewidth=2)
    ax.plot(iterations, z_analysis['max_trajectory'], 'r-', label='Max Z', linewidth=2)
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax.axhline(y=1, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Z Factor')
    ax.set_title('Quasiparticle Weights')
    ax.set_ylim(-0.1, 1.1)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Gap evolution
    gaps = [it['gap'] for it in analyzer.data['iterations']]
    ax = axes[1, 0]
    ax.plot(iterations, np.array(gaps) * 27.2114, 'g-', linewidth=2, marker='o', markersize=4)
    ax.set_xlabel('Iteration')
    ax.set_ylabel('HOMO-LUMO Gap (eV)')
    ax.set_title('Gap Evolution')
    ax.grid(True, alpha=0.3)

    # Memory usage
    memory = [it['memory_mb'] for it in analyzer.data['iterations']]
    ax = axes[1, 1]
    ax.plot(iterations, memory, 'r-', linewidth=2)
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Memory (MB)')
    ax.set_title('Memory Usage')
    ax.grid(True, alpha=0.3)

    if save_dir:
        save_dir = Path(save_dir)
        save_dir.mkdir(exist_ok=True)
        fig.savefig(save_dir / 'evgw_analysis.pdf', dpi=300, bbox_inches='tight')

    return fig


# Convenience factory functions
def create_monitor(config: Optional[Union[MonitorConfig, Dict[str, Any]]] = None) -> EvGWMonitor:
    """Create a new evGW monitor.

    Args:
        config: Configuration dict or MonitorConfig instance

    Returns:
        Configured EvGWMonitor instance
    """
    if isinstance(config, dict):
        config = MonitorConfig(**config)
    elif config is None:
        config = MonitorConfig()

    return EvGWMonitor(config)


@contextmanager
def monitor_context(config: Optional[Union[MonitorConfig, Dict[str, Any]]] = None):
    """Context manager for evGW monitoring.

    Args:
        config: Configuration dict or MonitorConfig instance

    Yields:
        EvGWMonitor instance
    """
    monitor = create_monitor(config)
    try:
        monitor.start()
        yield monitor
    except Exception as e:
        # Log the exception but ensure cleanup happens
        logger.error(f"Exception during monitoring: {e}")
        raise
    finally:
        # Always finalize, even if an exception occurred
        try:
            monitor.finalize()
        except Exception as e:
            logger.error(f"Error during monitor finalization: {e}")
            # Try emergency cleanup if finalize fails
            try:
                monitor._is_running = False
                if hasattr(monitor, '_plot_thread') and monitor._plot_thread:
                    monitor._update_queue.put(None)
            except:
                pass