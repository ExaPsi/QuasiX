"""
Logging configuration for QuasiX Python interface

Provides structured logging that integrates with the Rust core logging
and allows Python-level configuration and formatting.
"""

import json
import logging
import os
import sys
import time
from contextlib import contextmanager
from typing import Any, Dict, Optional


class QuasiXFormatter(logging.Formatter):
    """
    Custom formatter for QuasiX structured logging.
    
    Outputs logs in JSON format compatible with the Rust core logging.
    """
    
    def format(self, record: logging.LogRecord) -> str:
        """Format log record as JSON."""
        # Build base log entry
        log_entry = {
            "timestamp": self.formatTime(record, self.datefmt),
            "level": record.levelname,
            "target": f"quasix.{record.name}" if not record.name.startswith("quasix") else record.name,
            "message": record.getMessage(),
        }
        
        # Add optional fields if present
        if hasattr(record, "duration_ms"):
            log_entry["duration_ms"] = record.duration_ms
        
        if hasattr(record, "stage"):
            log_entry["stage"] = record.stage
            
        if hasattr(record, "context"):
            log_entry["context"] = record.context
            
        # Add thread info if in DEBUG mode
        if record.levelno <= logging.DEBUG:
            log_entry["thread_id"] = record.thread
            log_entry["thread_name"] = record.threadName
            
        # Add location info if available
        if hasattr(record, "pathname") and record.pathname:
            log_entry["file"] = f"{record.pathname}:{record.lineno}"
            log_entry["function"] = record.funcName
        
        return json.dumps(log_entry)


class QuasiXTextFormatter(logging.Formatter):
    """Human-readable text formatter for development."""
    
    def format(self, record: logging.LogRecord) -> str:
        """Format log record as human-readable text."""
        # Base format
        level = record.levelname[0]  # First letter only
        timestamp = self.formatTime(record, "%H:%M:%S.%f")[:-3]  # Millisecond precision
        
        message = f"[{timestamp}] {level} {record.name}: {record.getMessage()}"
        
        # Add timing info if present
        if hasattr(record, "duration_ms"):
            message += f" ({record.duration_ms:.1f}ms)"
            
        if hasattr(record, "stage"):
            message = f"[{timestamp}] {level} [{record.stage}] {record.name}: {record.getMessage()}"
            
        return message


def setup_logging(level: Optional[str] = None, format_type: Optional[str] = None) -> None:
    """
    Setup QuasiX logging configuration.
    
    Args:
        level: Log level (DEBUG, INFO, WARN, ERROR). If None, uses environment variable.
        format_type: Format type (json, text, pretty). If None, uses environment variable.
    """
    # Get configuration from environment or parameters
    log_level = level or os.environ.get("QUASIX_LOG_LEVEL", "INFO")
    log_format = format_type or os.environ.get("QUASIX_LOG_FORMAT", "json")
    
    # Convert to standard logging levels
    level_map = {
        "TRACE": logging.DEBUG,  # Python doesn't have TRACE
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARN": logging.WARNING,
        "WARNING": logging.WARNING,
        "ERROR": logging.ERROR,
    }
    
    numeric_level = level_map.get(log_level.upper(), logging.INFO)
    
    # Get or create logger
    logger = logging.getLogger("quasix")
    logger.setLevel(numeric_level)
    
    # Remove existing handlers
    logger.handlers.clear()
    
    # Create handler
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(numeric_level)
    
    # Set formatter based on format type
    if log_format.lower() == "json":
        formatter = QuasiXFormatter(datefmt="%Y-%m-%dT%H:%M:%S")
    elif log_format.lower() in ("text", "pretty"):
        formatter = QuasiXTextFormatter()
    else:
        # Default to JSON
        formatter = QuasiXFormatter(datefmt="%Y-%m-%dT%H:%M:%S")
    
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    
    # Also configure root logger to capture Rust logs if they come through
    if os.environ.get("QUASIX_CAPTURE_RUST_LOGS", "").lower() == "true":
        root_logger = logging.getLogger()
        root_logger.setLevel(numeric_level)
        root_handler = logging.StreamHandler(sys.stderr)
        root_handler.setFormatter(formatter)
        root_logger.addHandler(root_handler)
    
    # Log initialization
    logger.info("QuasiX Python logging initialized", extra={"context": {"level": log_level, "format": log_format}})


@contextmanager
def timed_stage(stage_name: str, logger: Optional[logging.Logger] = None):
    """
    Context manager for timing computational stages.
    
    Args:
        stage_name: Name of the stage being timed
        logger: Logger to use (defaults to quasix logger)
        
    Example:
        with timed_stage("frequency_integration"):
            # Expensive computation
            pass
    """
    if logger is None:
        logger = logging.getLogger("quasix")
    
    start_time = time.perf_counter()
    logger.info(f"Starting {stage_name}", extra={"stage": stage_name})
    
    try:
        yield
        duration_ms = (time.perf_counter() - start_time) * 1000
        logger.info(
            f"Completed {stage_name}",
            extra={"stage": stage_name, "duration_ms": duration_ms}
        )
    except Exception as e:
        duration_ms = (time.perf_counter() - start_time) * 1000
        logger.error(
            f"Failed {stage_name}: {e}",
            extra={"stage": stage_name, "duration_ms": duration_ms}
        )
        raise


class StageTimer:
    """Timer class for measuring stage durations."""
    
    def __init__(self, name: str, logger: Optional[logging.Logger] = None):
        """Initialize timer for a named stage."""
        self.name = name
        self.logger = logger or logging.getLogger("quasix")
        self.start_time = time.perf_counter()
        self.logger.info(f"Stage started: {name}", extra={"stage": name})
    
    def complete(self, success: bool = True, message: Optional[str] = None) -> float:
        """
        Complete the stage and log duration.
        
        Returns:
            Duration in milliseconds
        """
        duration_ms = (time.perf_counter() - self.start_time) * 1000
        
        if success:
            msg = message or f"Stage completed: {self.name}"
            self.logger.info(msg, extra={"stage": self.name, "duration_ms": duration_ms})
        else:
            msg = message or f"Stage failed: {self.name}"
            self.logger.error(msg, extra={"stage": self.name, "duration_ms": duration_ms})
        
        return duration_ms


def get_logger(name: str = "quasix") -> logging.Logger:
    """Get a QuasiX logger instance."""
    return logging.getLogger(name)


# Initialize logging on import if environment variables are set
if os.environ.get("QUASIX_LOG_LEVEL") or os.environ.get("QUASIX_LOG"):
    setup_logging()