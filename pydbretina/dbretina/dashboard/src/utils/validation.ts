/**
 * Frontend validation utilities for DBRetina dashboard.
 * Mirrors backend validation rules for immediate feedback.
 */

export interface ValidationResult {
  valid: boolean;
  error?: string;
}

/** Valid metrics matching backend VALID_METRICS */
export const VALID_METRICS = [
  "containment",
  "ochiai",
  "jaccard",
  "csi",
  "dice",
  "odds_ratio",
  "pvalue",
] as const;

export type MetricName = (typeof VALID_METRICS)[number];

/** Metric-specific cutoff ranges matching backend */
export const METRIC_CUTOFF_RANGES: Record<string, [number, number]> = {
  containment: [0, 100],
  ochiai: [0, 100],
  jaccard: [0, 100],
  csi: [0, 100],
  dice: [0, 100],
  odds_ratio: [0, Infinity],
  pvalue: [0, 1],
};

/**
 * Validate metric name
 */
export function validateMetric(metric: string): ValidationResult {
  if (!metric) {
    return { valid: false, error: "Metric is required" };
  }
  if (!VALID_METRICS.includes(metric as MetricName)) {
    return {
      valid: false,
      error: `Invalid metric "${metric}". Valid options: ${VALID_METRICS.join(", ")}`,
    };
  }
  return { valid: true };
}

/**
 * Validate cutoff value for a specific metric
 */
export function validateCutoff(
  cutoff: number | string,
  metric: string = "ochiai"
): ValidationResult {
  const numValue = typeof cutoff === "string" ? parseFloat(cutoff) : cutoff;

  if (isNaN(numValue)) {
    return { valid: false, error: "Cutoff must be a number" };
  }

  const [min, max] = METRIC_CUTOFF_RANGES[metric] || [0, 100];

  if (numValue < min) {
    return { valid: false, error: `Cutoff must be >= ${min}` };
  }

  if (max !== Infinity && numValue > max) {
    return { valid: false, error: `Cutoff must be <= ${max}` };
  }

  return { valid: true };
}

/**
 * Validate search query
 */
export function validateSearchQuery(query: string): ValidationResult {
  if (!query || query.trim().length === 0) {
    return { valid: false, error: "Search query is required" };
  }
  if (query.length < 1) {
    return { valid: false, error: "Search query must be at least 1 character" };
  }
  if (query.length > 200) {
    return { valid: false, error: "Search query must be less than 200 characters" };
  }
  return { valid: true };
}

/**
 * Validate SQL query for dangerous patterns
 */
const DANGEROUS_SQL_PATTERNS = [
  "DROP",
  "DELETE",
  "TRUNCATE",
  "ALTER",
  "CREATE",
  "INSERT",
  "UPDATE",
  "EXEC",
  "EXECUTE",
  "GRANT",
  "REVOKE",
  "--",
  "/*",
  "*/",
];

export function validateSQLSafety(query: string): ValidationResult {
  if (!query || query.trim().length === 0) {
    return { valid: false, error: "Query is required" };
  }

  const upper = query.toUpperCase();
  for (const pattern of DANGEROUS_SQL_PATTERNS) {
    if (upper.includes(pattern)) {
      return {
        valid: false,
        error: `Query contains blocked operation: ${pattern}`,
      };
    }
  }

  return { valid: true };
}

/**
 * Validate Cypher query for dangerous patterns
 */
const DANGEROUS_CYPHER_PATTERNS = [
  "DELETE",
  "DETACH DELETE",
  "CREATE",
  "MERGE",
  "SET",
  "REMOVE",
  "DROP",
  "CALL",
  "LOAD",
  "//",
];

export function validateCypherSafety(query: string): ValidationResult {
  if (!query || query.trim().length === 0) {
    return { valid: false, error: "Query is required" };
  }

  const upper = query.toUpperCase();
  for (const pattern of DANGEROUS_CYPHER_PATTERNS) {
    if (upper.includes(pattern)) {
      return {
        valid: false,
        error: `Query contains blocked operation: ${pattern}`,
      };
    }
  }

  return { valid: true };
}

/**
 * Validate pagination parameters
 */
export function validatePagination(
  limit: number,
  offset: number,
  maxLimit: number = 100000
): ValidationResult {
  if (limit < 1) {
    return { valid: false, error: "Limit must be at least 1" };
  }
  if (limit > maxLimit) {
    return { valid: false, error: `Limit must be at most ${maxLimit}` };
  }
  if (offset < 0) {
    return { valid: false, error: "Offset cannot be negative" };
  }
  return { valid: true };
}

/**
 * Validate hops parameter for neighborhood queries
 */
export function validateHops(hops: number): ValidationResult {
  if (!Number.isInteger(hops)) {
    return { valid: false, error: "Hops must be an integer" };
  }
  if (hops < 1) {
    return { valid: false, error: "Hops must be at least 1" };
  }
  if (hops > 5) {
    return { valid: false, error: "Hops must be at most 5 (to limit result size)" };
  }
  return { valid: true };
}

/**
 * Combined form validation helper
 */
export interface FormValidation {
  isValid: boolean;
  errors: Record<string, string>;
}

export function validateForm(
  fields: Record<string, { value: unknown; validator: (v: unknown) => ValidationResult }>
): FormValidation {
  const errors: Record<string, string> = {};
  let isValid = true;

  for (const [fieldName, { value, validator }] of Object.entries(fields)) {
    const result = validator(value);
    if (!result.valid) {
      errors[fieldName] = result.error || "Invalid value";
      isValid = false;
    }
  }

  return { isValid, errors };
}
