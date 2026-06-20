"""Custom exceptions for DBRetina REST API.

This module provides a structured exception hierarchy for the REST API,
enabling consistent error responses with appropriate HTTP status codes
and user-friendly messages.
"""

import re
from typing import Any, Optional


class DBRetinaAPIError(Exception):
    """Base exception for all DBRetina API errors.

    Attributes:
        status_code: HTTP status code for the error response.
        detail: User-friendly error message.
        error_code: Machine-readable error code for client handling.
        context: Additional context data for debugging.
    """

    status_code: int = 500
    error_code: str = "INTERNAL_ERROR"

    def __init__(
        self,
        detail: str = "An internal error occurred",
        context: Optional[dict[str, Any]] = None,
    ):
        self.detail = detail
        self.context = context or {}
        super().__init__(detail)

    def to_dict(self) -> dict[str, Any]:
        """Convert exception to API response format."""
        response = {
            "error": True,
            "error_code": self.error_code,
            "detail": self.detail,
        }
        if self.context:
            response["context"] = self.context
        return response


class DataNotFoundError(DBRetinaAPIError):
    """Raised when requested data (group, pair, etc.) is not found."""

    status_code = 404
    error_code = "DATA_NOT_FOUND"

    def __init__(
        self,
        detail: str = "The requested data was not found",
        resource_type: Optional[str] = None,
        resource_id: Optional[str] = None,
    ):
        context = {}
        if resource_type:
            context["resource_type"] = resource_type
        if resource_id:
            context["resource_id"] = resource_id
        super().__init__(detail, context)


class ValidationError(DBRetinaAPIError):
    """Raised when input validation fails."""

    status_code = 400
    error_code = "VALIDATION_ERROR"

    def __init__(
        self,
        detail: str = "Input validation failed",
        field: Optional[str] = None,
        value: Optional[Any] = None,
        allowed_values: Optional[list] = None,
    ):
        context = {}
        if field:
            context["field"] = field
        if value is not None:
            context["value"] = str(value)
        if allowed_values:
            context["allowed_values"] = allowed_values
        super().__init__(detail, context)


class QueryTimeoutError(DBRetinaAPIError):
    """Raised when a query exceeds the timeout limit."""

    status_code = 408
    error_code = "QUERY_TIMEOUT"

    def __init__(
        self,
        detail: str = "Query timed out",
        timeout_seconds: Optional[float] = None,
        query_type: Optional[str] = None,
    ):
        context = {}
        if timeout_seconds is not None:
            context["timeout_seconds"] = timeout_seconds
        if query_type:
            context["query_type"] = query_type
        super().__init__(detail, context)


class DataTooLargeError(DBRetinaAPIError):
    """Raised when requested data exceeds size limits."""

    status_code = 413
    error_code = "DATA_TOO_LARGE"

    def __init__(
        self,
        detail: str = "Requested data exceeds size limit",
        requested_size: Optional[int] = None,
        max_size: Optional[int] = None,
        suggestion: Optional[str] = None,
    ):
        context = {}
        if requested_size is not None:
            context["requested_size"] = requested_size
        if max_size is not None:
            context["max_size"] = max_size
        if suggestion:
            context["suggestion"] = suggestion
        super().__init__(detail, context)


class QuerySyntaxError(DBRetinaAPIError):
    """Raised when a SQL query has syntax errors."""

    status_code = 400
    error_code = "QUERY_SYNTAX_ERROR"

    def __init__(
        self,
        detail: str = "Query syntax error",
        query_type: Optional[str] = None,
        query: Optional[str] = None,
        position: Optional[int] = None,
    ):
        context = {}
        if query_type:
            context["query_type"] = query_type
        if query:
            # Truncate long queries
            context["query"] = query[:200] + "..." if len(query) > 200 else query
        if position is not None:
            context["position"] = position
        super().__init__(detail, context)


class UnsafeQueryError(DBRetinaAPIError):
    """Raised when a query contains potentially dangerous operations."""

    status_code = 400
    error_code = "UNSAFE_QUERY"

    def __init__(
        self,
        detail: str = "Query contains unsafe operations",
        operation: Optional[str] = None,
    ):
        context = {}
        if operation:
            context["blocked_operation"] = operation
        super().__init__(detail, context)


class FeatureNotAvailableError(DBRetinaAPIError):
    """Raised when a feature requires additional setup (e.g., gene index)."""

    status_code = 400
    error_code = "FEATURE_NOT_AVAILABLE"

    def __init__(
        self,
        detail: str = "This feature is not available",
        feature: Optional[str] = None,
        requirement: Optional[str] = None,
    ):
        context = {}
        if feature:
            context["feature"] = feature
        if requirement:
            context["requirement"] = requirement
        super().__init__(detail, context)


class AlgorithmError(DBRetinaAPIError):
    """Raised when a graph algorithm fails."""

    status_code = 500
    error_code = "ALGORITHM_ERROR"

    def __init__(
        self,
        detail: str = "Algorithm execution failed",
        algorithm: Optional[str] = None,
        reason: Optional[str] = None,
    ):
        context = {}
        if algorithm:
            context["algorithm"] = algorithm
        if reason:
            context["reason"] = reason
        super().__init__(detail, context)


# ── Utility functions ─────────────────────────────────────────────

DANGEROUS_SQL_PATTERNS = [
    "DROP", "DELETE", "TRUNCATE", "ALTER", "CREATE", "INSERT",
    "UPDATE", "EXEC", "EXECUTE", "GRANT", "REVOKE",
]


def _strip_string_literals(query: str) -> str:
    """Remove string literal contents so they don't trigger keyword checks."""
    return re.sub(r"'[^']*'|\"[^\"]*\"", "''", query)


def validate_sql_safety(query: str) -> None:
    """Check if SQL query is safe for read-only execution.

    Strips string literals and uses word-boundary matching so values
    like 'updated_pathway' don't trigger the UPDATE check.

    Raises:
        UnsafeQueryError: If query contains dangerous patterns.
    """
    stripped = _strip_string_literals(query)
    for pattern in DANGEROUS_SQL_PATTERNS:
        if re.search(r"\b" + re.escape(pattern) + r"\b", stripped, re.IGNORECASE):
            raise UnsafeQueryError(
                detail=f"SQL query contains blocked operation: {pattern}",
                operation=pattern,
            )


def validate_metric(metric: str, valid_metrics: tuple[str, ...]) -> None:
    """Validate that a metric name is valid.

    Raises:
        ValidationError: If metric is not in valid_metrics.
    """
    if metric not in valid_metrics:
        raise ValidationError(
            detail=f"Unknown metric '{metric}'",
            field="metric",
            value=metric,
            allowed_values=list(valid_metrics),
        )


def validate_cutoff(cutoff: float, min_val: float = 0.0, max_val: float = 100.0) -> None:
    """Validate that a cutoff value is within range.

    Raises:
        ValidationError: If cutoff is out of range.
    """
    if cutoff < min_val or cutoff > max_val:
        raise ValidationError(
            detail=f"Cutoff must be between {min_val} and {max_val}",
            field="cutoff",
            value=cutoff,
            allowed_values=[f"{min_val}-{max_val}"],
        )
