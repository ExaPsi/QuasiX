"""
Basic tests for QuasiX Python bindings - S1-2 verification
"""

import pytest
import sys
import threading
import time

def test_import():
    """Test that quasix can be imported"""
    import quasix
    assert quasix is not None
    
def test_version():
    """Test version function returns string"""
    import quasix
    version = quasix.version()
    assert isinstance(version, str)
    assert version == "0.6.0"
    
def test_metadata():
    """Test metadata function returns dict with expected keys"""
    import quasix
    meta = quasix.metadata()
    assert isinstance(meta, dict)
    assert "version" in meta
    assert "name" in meta
    assert "authors" in meta
    assert "description" in meta
    assert "rust_version" in meta
    assert meta["version"] == "0.6.0"
    assert meta["name"] == "quasix"
    
def test_noop_kernel():
    """Test no-op kernel executes successfully"""
    import quasix
    result = quasix.noop_kernel()
    assert isinstance(result, dict)
    assert result["status"] == "success"
    assert "message" in result
    assert result["version"] == "0.6.0"
    
def test_module_attributes():
    """Test module has expected attributes"""
    import quasix
    attrs = dir(quasix)
    assert "version" in attrs
    assert "metadata" in attrs
    assert "noop_kernel" in attrs
    assert "__version__" in attrs
    
def test_version_consistency():
    """Test that __version__ matches version()"""
    import quasix
    assert quasix.__version__ == quasix.version()
    
def test_gil_release_threading():
    """Test that GIL is released during calls (basic threading test)"""
    import quasix
    
    results = []
    
    def worker():
        result = quasix.noop_kernel()
        results.append(result)
        return result
    
    # Test parallel execution
    start = time.time()
    threads = [threading.Thread(target=worker) for _ in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    parallel_time = time.time() - start
    
    # Clear results
    results.clear()
    
    # Test sequential execution
    start = time.time()
    for _ in range(4):
        worker()
    sequential_time = time.time() - start
    
    # All results should be successful
    assert len(results) == 4
    for result in results:
        assert result["status"] == "success"
    
    # Note: The time comparison might not show difference for no-op
    # but we're verifying that threading works without deadlocks
    print(f"Parallel: {parallel_time:.3f}s, Sequential: {sequential_time:.3f}s")
    
def test_import_rust_extension():
    """Test that the Rust extension module can be imported directly"""
    import quasix.quasix as rust_module
    assert rust_module is not None
    assert hasattr(rust_module, "version")
    assert hasattr(rust_module, "metadata")
    assert hasattr(rust_module, "noop_kernel")
    
if __name__ == "__main__":
    # Run basic smoke tests
    print("Running QuasiX S1-2 verification tests...")
    test_import()
    print("✓ Import test passed")
    test_version()
    print("✓ Version test passed")
    test_metadata()
    print("✓ Metadata test passed")
    test_noop_kernel()
    print("✓ No-op kernel test passed")
    test_module_attributes()
    print("✓ Module attributes test passed")
    test_version_consistency()
    print("✓ Version consistency test passed")
    test_gil_release_threading()
    print("✓ GIL release threading test passed")
    test_import_rust_extension()
    print("✓ Rust extension import test passed")
    print("\nAll S1-2 verification tests passed!")