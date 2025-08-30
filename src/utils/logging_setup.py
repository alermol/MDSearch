"""Logging configuration utilities."""

import logging
import time
import json
from typing import Optional

__all__ = ["JsonFormatter", "setup_logger"]


class JsonFormatter(logging.Formatter):
    """Custom JSON formatter for structured logging."""

    def format(self, record: logging.LogRecord) -> str:
        """Format log record as JSON string."""
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
    """Configure and return a logger with specified settings."""
    logger = logging.getLogger(name)

    if not logger.handlers:
        handler = logging.StreamHandler()

        if format_type == "json":
            handler.setFormatter(JsonFormatter())
        else:
            handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))

        logger.addHandler(handler)
        logger.propagate = False

    level_name = (level or ("INFO" if verbose else "WARNING")).upper()
    log_level = getattr(logging, level_name, logging.INFO)
    logger.setLevel(log_level)

    return logger
