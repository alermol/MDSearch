"""Logging configuration utilities."""

import logging
import time
import json
from typing import Optional

__all__ = ["JsonFormatter", "setup_logger"]


class JsonFormatter(logging.Formatter):
    """Custom JSON formatter for structured logging.

    Formats log records as JSON with timestamp, level, and message.

    Example:
        >>> formatter = JsonFormatter()
        >>> record = logging.LogRecord(
        ...     name="test", level=logging.INFO, pathname="", lineno=0,
        ...     msg="Test message", args=(), exc_info=None
        ... )
        >>> record.created = time.time()
        >>> json_output = formatter.format(record)
        >>> print(json_output)
        {"level": "INFO", "message": "Test message", "time": "2024-01-01T12:00:00Z"}
    """

    def format(self, record: logging.LogRecord) -> str:
        """Format log record as JSON string.

        Args:
            record: LogRecord to format

        Returns:
            JSON-formatted log string

        Example:
            >>> formatter = JsonFormatter()
            >>> record = logging.LogRecord(
            ...     name="test", level=logging.WARNING, pathname="", lineno=0,
            ...     msg="Warning message", args=(), exc_info=None
            ... )
            >>> record.created = time.time()
            >>> json_output = formatter.format(record)
            >>> print(json_output)
            {"level": "WARNING", "message": "Warning message", "time": "2024-01-01T12:00:00Z"}
        """
        payload = {
            "level": record.levelname,
            "message": record.getMessage(),
            "time": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(record.created)),
        }
        return json.dumps(payload, ensure_ascii=False)


def setup_logger(
    name: str,
    level: Optional[str] = None,
    format_type: str = "text",
    verbose: bool = True,
) -> logging.Logger:
    """Configure and return a logger with specified settings.

    Args:
        name: Logger name
        level: Logging level (DEBUG, INFO, WARNING, ERROR)
        format_type: Format type ("text" or "json")
        verbose: Whether to enable verbose logging

    Returns:
        Configured logger instance

    Example:
        >>> # Create text logger with INFO level
        >>> logger = setup_logger("my_app", level="INFO", format_type="text")
        >>> logger.info("Application started")
        [INFO] Application started

        >>> # Create JSON logger with DEBUG level
        >>> json_logger = setup_logger("api", level="DEBUG", format_type="json")
        >>> json_logger.debug("API request received")
        {"level": "DEBUG", "message": "API request received", "time": "2024-01-01T12:00:00Z"}

        >>> # Create quiet logger (WARNING level)
        >>> quiet_logger = setup_logger("background", verbose=False)
        >>> quiet_logger.info("This won't be logged")
        >>> quiet_logger.warning("This will be logged")
        [WARNING] This will be logged
    """
    logger = logging.getLogger(name)

    # Avoid duplicate handlers if multiple instances
    if not logger.handlers:
        handler = logging.StreamHandler()

        if format_type == "json":
            handler.setFormatter(JsonFormatter())
        else:
            handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))

        logger.addHandler(handler)
        logger.propagate = False

    # Determine level
    level_name = (level or ("INFO" if verbose else "WARNING")).upper()
    log_level = getattr(logging, level_name, logging.INFO)
    logger.setLevel(log_level)

    return logger
