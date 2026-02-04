import { useState, useCallback, useMemo } from "react";
import { ValidationResult } from "../utils/validation";

export interface FieldState {
  value: string;
  error: string | null;
  touched: boolean;
  dirty: boolean;
}

export interface UseValidationOptions<T extends Record<string, string>> {
  initialValues: T;
  validators: Partial<Record<keyof T, (value: string) => ValidationResult>>;
  validateOnChange?: boolean;
  validateOnBlur?: boolean;
}

export interface UseValidationReturn<T extends Record<string, string>> {
  values: T;
  errors: Partial<Record<keyof T, string>>;
  touched: Partial<Record<keyof T, boolean>>;
  isValid: boolean;
  isDirty: boolean;
  setValue: (field: keyof T, value: string) => void;
  setTouched: (field: keyof T) => void;
  validateField: (field: keyof T) => boolean;
  validateAll: () => boolean;
  reset: () => void;
  getFieldProps: (field: keyof T) => {
    value: string;
    onChange: (e: React.ChangeEvent<HTMLInputElement | HTMLSelectElement | HTMLTextAreaElement>) => void;
    onBlur: () => void;
    error?: string;
  };
}

/**
 * Hook for managing form validation state.
 * Provides field-level and form-level validation with dirty/touched tracking.
 */
export function useValidation<T extends Record<string, string>>(
  options: UseValidationOptions<T>
): UseValidationReturn<T> {
  const { initialValues, validators, validateOnChange = true, validateOnBlur = true } = options;

  const [values, setValues] = useState<T>(initialValues);
  const [errors, setErrors] = useState<Partial<Record<keyof T, string>>>({});
  const [touched, setTouchedFields] = useState<Partial<Record<keyof T, boolean>>>({});

  const validateField = useCallback(
    (field: keyof T): boolean => {
      const validator = validators[field];
      if (!validator) return true;

      const result = validator(values[field]);
      setErrors((prev) => ({
        ...prev,
        [field]: result.valid ? undefined : result.error,
      }));
      return result.valid;
    },
    [validators, values]
  );

  const validateAll = useCallback((): boolean => {
    let isValid = true;
    const newErrors: Partial<Record<keyof T, string>> = {};

    for (const field of Object.keys(values) as Array<keyof T>) {
      const validator = validators[field];
      if (validator) {
        const result = validator(values[field]);
        if (!result.valid) {
          newErrors[field] = result.error;
          isValid = false;
        }
      }
    }

    setErrors(newErrors);
    // Mark all fields as touched
    const allTouched = Object.keys(values).reduce(
      (acc, key) => ({ ...acc, [key]: true }),
      {} as Partial<Record<keyof T, boolean>>
    );
    setTouchedFields(allTouched);

    return isValid;
  }, [validators, values]);

  const setValue = useCallback(
    (field: keyof T, value: string) => {
      setValues((prev) => ({ ...prev, [field]: value }));
      if (validateOnChange) {
        const validator = validators[field];
        if (validator) {
          const result = validator(value);
          setErrors((prev) => ({
            ...prev,
            [field]: result.valid ? undefined : result.error,
          }));
        }
      }
    },
    [validateOnChange, validators]
  );

  const setTouched = useCallback(
    (field: keyof T) => {
      setTouchedFields((prev) => ({ ...prev, [field]: true }));
      if (validateOnBlur) {
        validateField(field);
      }
    },
    [validateOnBlur, validateField]
  );

  const reset = useCallback(() => {
    setValues(initialValues);
    setErrors({});
    setTouchedFields({});
  }, [initialValues]);

  const isValid = useMemo(() => {
    return Object.keys(errors).every(
      (key) => errors[key as keyof T] === undefined
    );
  }, [errors]);

  const isDirty = useMemo(() => {
    return Object.keys(values).some(
      (key) => values[key as keyof T] !== initialValues[key as keyof T]
    );
  }, [values, initialValues]);

  const getFieldProps = useCallback(
    (field: keyof T) => ({
      value: values[field],
      onChange: (
        e: React.ChangeEvent<HTMLInputElement | HTMLSelectElement | HTMLTextAreaElement>
      ) => setValue(field, e.target.value),
      onBlur: () => setTouched(field),
      error: touched[field] ? errors[field] : undefined,
    }),
    [values, errors, touched, setValue, setTouched]
  );

  return {
    values,
    errors,
    touched,
    isValid,
    isDirty,
    setValue,
    setTouched,
    validateField,
    validateAll,
    reset,
    getFieldProps,
  };
}

/**
 * Simple single-field validation hook
 */
export function useFieldValidation(
  initialValue: string,
  validator: (value: string) => ValidationResult
) {
  const [value, setValueState] = useState(initialValue);
  const [error, setError] = useState<string | null>(null);
  const [touched, setTouched] = useState(false);

  const setValue = useCallback(
    (newValue: string) => {
      setValueState(newValue);
      const result = validator(newValue);
      setError(result.valid ? null : result.error || "Invalid");
    },
    [validator]
  );

  const blur = useCallback(() => {
    setTouched(true);
    const result = validator(value);
    setError(result.valid ? null : result.error || "Invalid");
  }, [validator, value]);

  const reset = useCallback(() => {
    setValueState(initialValue);
    setError(null);
    setTouched(false);
  }, [initialValue]);

  return {
    value,
    error: touched ? error : null,
    isValid: error === null,
    setValue,
    blur,
    reset,
    getInputProps: () => ({
      value,
      onChange: (e: React.ChangeEvent<HTMLInputElement>) => setValue(e.target.value),
      onBlur: blur,
    }),
  };
}
