/**
 * Custom React hooks for DBRetina dashboard.
 */

export {
  useDebounce,
  useDebouncedCallback,
  useDebounceWithFlush,
} from "./useDebounce";

export {
  useValidation,
  useFieldValidation,
  type FieldState,
  type UseValidationOptions,
  type UseValidationReturn,
} from "./useValidation";

export {
  useQueryHistory,
  SQL_TEMPLATES,
  CYPHER_TEMPLATES,
  type QueryHistoryItem,
} from "./useQueryHistory";
