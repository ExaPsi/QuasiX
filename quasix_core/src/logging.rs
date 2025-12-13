//! Logging and tracing infrastructure for QuasiX
//!
//! Provides structured logging with JSON output, configurable levels,
//! and stage timing functionality.

use std::time::Instant;
use tracing::info;
use tracing_subscriber::{
    fmt::{self, format::FmtSpan},
    layer::SubscriberExt,
    util::SubscriberInitExt,
    EnvFilter,
};

/// Initialize the QuasiX logger with environment-based configuration
///
/// # Environment Variables
/// - `QUASIX_LOG`: Log level (TRACE, DEBUG, INFO, WARN, ERROR)
/// - `QUASIX_LOG_FORMAT`: Output format (json, pretty, compact)
/// - `RUST_LOG`: Standard Rust log level (fallback)
pub fn init_logger() -> anyhow::Result<()> {
    // Check if already initialized
    static INIT: std::sync::Once = std::sync::Once::new();
    let mut result = Ok(());

    INIT.call_once(|| {
        // Create environment filter from QUASIX_LOG or RUST_LOG
        let env_filter = std::env::var("QUASIX_LOG")
            .or_else(|_| std::env::var("RUST_LOG"))
            .unwrap_or_else(|_| "quasix_core=info".to_string());

        let filter =
            EnvFilter::try_new(&env_filter).unwrap_or_else(|_| EnvFilter::new("quasix_core=info"));

        // Check format preference
        let format = std::env::var("QUASIX_LOG_FORMAT").unwrap_or_else(|_| "json".to_string());

        match format.as_str() {
            "json" => {
                // JSON formatted output
                let fmt_layer = fmt::layer()
                    .json()
                    .with_current_span(true)
                    .with_span_list(false)
                    .with_target(true)
                    .with_level(true)
                    .with_thread_ids(false)
                    .with_thread_names(false)
                    .with_file(false)
                    .with_line_number(false)
                    .with_span_events(FmtSpan::CLOSE);

                if let Err(e) = tracing_subscriber::registry()
                    .with(filter)
                    .with(fmt_layer)
                    .try_init()
                {
                    result = Err(anyhow::anyhow!("Failed to initialize logger: {}", e));
                }
            }
            "pretty" => {
                // Pretty formatted output for development
                let fmt_layer = fmt::layer()
                    .pretty()
                    .with_target(true)
                    .with_level(true)
                    .with_thread_ids(false)
                    .with_thread_names(false)
                    .with_span_events(FmtSpan::CLOSE);

                if let Err(e) = tracing_subscriber::registry()
                    .with(filter)
                    .with(fmt_layer)
                    .try_init()
                {
                    result = Err(anyhow::anyhow!("Failed to initialize logger: {}", e));
                }
            }
            _ => {
                // Compact format (default fallback)
                let fmt_layer = fmt::layer()
                    .compact()
                    .with_target(true)
                    .with_level(true)
                    .with_thread_ids(false)
                    .with_thread_names(false)
                    .with_span_events(FmtSpan::CLOSE);

                if let Err(e) = tracing_subscriber::registry()
                    .with(filter)
                    .with(fmt_layer)
                    .try_init()
                {
                    result = Err(anyhow::anyhow!("Failed to initialize logger: {}", e));
                }
            }
        }

        if result.is_ok() {
            info!("QuasiX logger initialized with format: {}", format);
        }
    });

    result
}

/// Timer for measuring stage durations
pub struct StageTimer {
    name: String,
    start: Instant,
}

impl StageTimer {
    /// Create a new stage timer
    #[must_use]
    pub fn new(name: impl Into<String>) -> Self {
        let name = name.into();
        info!(stage = %name, "Stage started");
        Self {
            name,
            start: Instant::now(),
        }
    }

    /// Complete the stage and log the duration
    pub fn complete(self) {
        let duration_ms = self.start.elapsed().as_secs_f64() * 1000.0;
        info!(
            stage = %self.name,
            duration_ms = duration_ms,
            "Stage completed"
        );
    }

    /// Complete with a result
    pub fn complete_with_result<T>(self, result: &Result<T, impl std::fmt::Display>) {
        let duration_ms = self.start.elapsed().as_secs_f64() * 1000.0;
        match result {
            Ok(_) => {
                info!(
                    stage = %self.name,
                    duration_ms = duration_ms,
                    status = "success",
                    "Stage completed successfully"
                );
            }
            Err(e) => {
                tracing::error!(
                    stage = %self.name,
                    duration_ms = duration_ms,
                    status = "error",
                    error = %e,
                    "Stage failed"
                );
            }
        }
    }
}

/// Macro for timing a code block
#[macro_export]
macro_rules! timed_stage {
    ($stage:expr, $block:block) => {{
        let _timer = $crate::logging::StageTimer::new($stage);
        let result = $block;
        _timer.complete();
        result
    }};
}

/// Macro for timing a code block with result handling
#[macro_export]
macro_rules! timed_stage_result {
    ($stage:expr, $block:block) => {{
        let _timer = $crate::logging::StageTimer::new($stage);
        let result = $block;
        _timer.complete_with_result(&result);
        result
    }};
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_logger_initialization() {
        // Set test environment
        std::env::set_var("QUASIX_LOG", "debug");
        std::env::set_var("QUASIX_LOG_FORMAT", "json");

        // Initialize logger
        let result = init_logger();
        assert!(result.is_ok());

        // Test logging
        info!("Test log message");
        tracing::debug!("Debug message");
        tracing::warn!("Warning message");
    }

    #[test]
    fn test_stage_timer() {
        init_logger().ok();

        let timer = StageTimer::new("test_stage");
        std::thread::sleep(std::time::Duration::from_millis(10));
        timer.complete();
    }

    #[test]
    fn test_timed_stage_macro() {
        init_logger().ok();

        let result = timed_stage!("macro_test", {
            std::thread::sleep(std::time::Duration::from_millis(10));
            42
        });

        assert_eq!(result, 42);
    }
}
