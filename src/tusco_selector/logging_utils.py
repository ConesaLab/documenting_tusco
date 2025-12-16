"""Centralized logging utilities for TUSCO selector."""

import logging
import sys
from pathlib import Path
from typing import Optional

# Custom log levels
TRACE = 5
VERBOSE = 15

# Register custom levels
logging.addLevelName(TRACE, "TRACE")
logging.addLevelName(VERBOSE, "VERBOSE")


class ColoredFormatter(logging.Formatter):
    """Colored formatter for console output."""
    
    # ANSI color codes
    COLORS = {
        'TRACE': '\033[90m',      # Gray
        'DEBUG': '\033[36m',       # Cyan
        'VERBOSE': '\033[34m',     # Blue
        'INFO': '\033[32m',        # Green
        'WARNING': '\033[33m',     # Yellow
        'ERROR': '\033[31m',       # Red
        'CRITICAL': '\033[35m',    # Magenta
        'RESET': '\033[0m',
    }
    
    def __init__(self, fmt: Optional[str] = None, use_colors: bool = True):
        super().__init__(fmt)
        self.use_colors = use_colors and sys.stderr.isatty()
    
    def format(self, record: logging.LogRecord) -> str:
        if self.use_colors:
            levelname = record.levelname
            if levelname in self.COLORS:
                record.levelname = f"{self.COLORS[levelname]}{levelname}{self.COLORS['RESET']}"
        return super().format(record)


def init_logging(
    level: int = logging.INFO,
    log_file: Optional[Path] = None,
    use_colors: bool = True,
    format_string: Optional[str] = None
) -> None:
    """Initialize logging configuration.
    
    Args:
        level: Logging level (can use TRACE or VERBOSE custom levels)
        log_file: Optional file to write logs to
        use_colors: Whether to use colored output for console
        format_string: Custom format string (defaults to sensible format)
    """
    if format_string is None:
        format_string = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    
    # Clear existing handlers
    root_logger.handlers.clear()
    
    # Console handler with colored formatter
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(level)
    console_formatter = ColoredFormatter(format_string, use_colors=use_colors)
    console_handler.setFormatter(console_formatter)
    root_logger.addHandler(console_handler)
    
    # File handler if specified
    if log_file is not None:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_formatter = logging.Formatter(format_string)
        file_handler.setFormatter(file_formatter)
        root_logger.addHandler(file_handler)


def get_logger(name: str) -> logging.Logger:
    """Get a logger instance with custom methods.
    
    Args:
        name: Logger name (typically __name__)
        
    Returns:
        Logger instance with trace and verbose methods
    """
    logger = logging.getLogger(name)
    
    # Add custom level methods if not already present
    if not hasattr(logger, 'trace'):
        def trace(msg, *args, **kwargs):
            if logger.isEnabledFor(TRACE):
                logger._log(TRACE, msg, args, **kwargs)
        logger.trace = trace
    
    if not hasattr(logger, 'verbose'):
        def verbose(msg, *args, **kwargs):
            if logger.isEnabledFor(VERBOSE):
                logger._log(VERBOSE, msg, args, **kwargs)
        logger.verbose = verbose
    
    return logger


def set_log_level(level: int) -> None:
    """Set the logging level for all handlers.
    
    Args:
        level: New logging level
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    for handler in root_logger.handlers:
        handler.setLevel(level)


# Convenience level setting functions
def set_trace() -> None:
    """Set logging to TRACE level."""
    set_log_level(TRACE)


def set_verbose() -> None:
    """Set logging to VERBOSE level."""
    set_log_level(VERBOSE)


def set_debug() -> None:
    """Set logging to DEBUG level."""
    set_log_level(logging.DEBUG)


def set_info() -> None:
    """Set logging to INFO level."""
    set_log_level(logging.INFO)


__all__ = [
    'TRACE',
    'VERBOSE',
    'init_logging',
    'get_logger',
    'set_log_level',
    'set_trace',
    'set_verbose',
    'set_debug',
    'set_info',
    'ColoredFormatter',
]