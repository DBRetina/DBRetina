import { useState, useCallback } from "react";
import { useDashboard } from "../../../state/context";
import { VALID_METRICS } from "../../../utils/validation";
import { genId } from "../../../utils/uid";

export interface FilterCondition {
  id: string;
  metric: string;
  operator: string;
  value: number | [number, number];
}

interface FilterConditionRowProps {
  condition: FilterCondition;
  onChange: (condition: FilterCondition) => void;
  onRemove: () => void;
  canRemove: boolean;
}

function FilterConditionRow({ condition, onChange, onRemove, canRemove }: FilterConditionRowProps) {
  const operators = [
    { value: ">=", label: ">=" },
    { value: "<=", label: "<=" },
    { value: ">", label: ">" },
    { value: "<", label: "<" },
    { value: "==", label: "=" },
    { value: "!=", label: "!=" },
    { value: "between", label: "between" },
  ];

  const isBetween = condition.operator === "between";
  const betweenValue = Array.isArray(condition.value) ? condition.value : [0, 100];

  return (
    <div
      style={{
        display: "flex",
        flexDirection: "column",
        gap: 8,
        padding: 10,
        background: "var(--bg-primary)",
        borderRadius: 4,
        border: "1px solid var(--border-color)",
      }}
    >
      <div style={{ display: "flex", gap: 8, alignItems: "center" }}>
        {/* Metric selector */}
        <select
          value={condition.metric}
          onChange={(e) => onChange({ ...condition, metric: e.target.value })}
          style={{
            flex: 1,
            padding: "4px 8px",
            border: "1px solid var(--border-color)",
            borderRadius: 4,
            fontSize: 12,
            background: "white",
          }}
        >
          {VALID_METRICS.map((metric) => (
            <option key={metric} value={metric}>
              {metric}
            </option>
          ))}
        </select>

        {/* Operator selector */}
        <select
          value={condition.operator}
          onChange={(e) => {
            const newOp = e.target.value;
            const newValue: number | [number, number] = newOp === "between" ? [0, 100] as [number, number] : (Array.isArray(condition.value) ? condition.value[0] : condition.value);
            onChange({ ...condition, operator: newOp, value: newValue });
          }}
          style={{
            width: 80,
            padding: "4px 8px",
            border: "1px solid var(--border-color)",
            borderRadius: 4,
            fontSize: 12,
            background: "white",
          }}
        >
          {operators.map((op) => (
            <option key={op.value} value={op.value}>
              {op.label}
            </option>
          ))}
        </select>

        {/* Remove button */}
        {canRemove && (
          <button
            onClick={onRemove}
            style={{
              padding: "4px 8px",
              border: "none",
              background: "transparent",
              color: "var(--danger)",
              cursor: "pointer",
              fontSize: 14,
            }}
            title="Remove condition"
          >
            &times;
          </button>
        )}
      </div>

      {/* Value input */}
      {isBetween ? (
        <div style={{ display: "flex", gap: 8, alignItems: "center" }}>
          <input
            type="number"
            value={betweenValue[0]}
            onChange={(e) =>
              onChange({
                ...condition,
                value: [parseFloat(e.target.value) || 0, betweenValue[1]],
              })
            }
            style={{
              flex: 1,
              padding: "4px 8px",
              border: "1px solid var(--border-color)",
              borderRadius: 4,
              fontSize: 12,
            }}
            placeholder="Min"
          />
          <span style={{ fontSize: 11, color: "var(--text-secondary)" }}>to</span>
          <input
            type="number"
            value={betweenValue[1]}
            onChange={(e) =>
              onChange({
                ...condition,
                value: [betweenValue[0], parseFloat(e.target.value) || 0],
              })
            }
            style={{
              flex: 1,
              padding: "4px 8px",
              border: "1px solid var(--border-color)",
              borderRadius: 4,
              fontSize: 12,
            }}
            placeholder="Max"
          />
        </div>
      ) : (
        <input
          type="number"
          value={Array.isArray(condition.value) ? condition.value[0] : condition.value}
          onChange={(e) =>
            onChange({ ...condition, value: parseFloat(e.target.value) || 0 })
          }
          style={{
            width: "100%",
            padding: "4px 8px",
            border: "1px solid var(--border-color)",
            borderRadius: 4,
            fontSize: 12,
          }}
          placeholder="Value"
        />
      )}
    </div>
  );
}

interface FilterBuilderProps {
  onApply: (filters: FilterCondition[], logic: "AND" | "OR") => void;
  onClear: () => void;
}

export default function FilterBuilder({ onApply, onClear }: FilterBuilderProps) {
  const { state } = useDashboard();
  const [conditions, setConditions] = useState<FilterCondition[]>([
    { id: genId(), metric: "ochiai", operator: ">=", value: 50 },
  ]);
  const [logic, setLogic] = useState<"AND" | "OR">("AND");
  const [expanded, setExpanded] = useState(false);

  // Presets
  const presets = [
    {
      name: "High Confidence",
      conditions: [
        { metric: "ochiai", operator: ">=", value: 70 },
        { metric: "jaccard", operator: ">=", value: 50 },
      ],
      logic: "AND" as const,
    },
    {
      name: "Significant Only",
      conditions: [{ metric: "pvalue", operator: "<=", value: 0.05 }],
      logic: "AND" as const,
    },
    {
      name: "Strong Overlap",
      conditions: [{ metric: "containment", operator: ">=", value: 80 }],
      logic: "AND" as const,
    },
  ];

  const addCondition = useCallback(() => {
    setConditions((prev) => [
      ...prev,
      { id: genId(), metric: "ochiai", operator: ">=", value: 50 },
    ]);
  }, []);

  const updateCondition = useCallback((id: string, updated: FilterCondition) => {
    setConditions((prev) => prev.map((c) => (c.id === id ? updated : c)));
  }, []);

  const removeCondition = useCallback((id: string) => {
    setConditions((prev) => prev.filter((c) => c.id !== id));
  }, []);

  const applyPreset = useCallback((preset: (typeof presets)[0]) => {
    setConditions(
      preset.conditions.map((c) => ({
        id: genId(),
        metric: c.metric,
        operator: c.operator,
        value: c.value,
      }))
    );
    setLogic(preset.logic);
  }, []);

  const handleApply = useCallback(() => {
    onApply(conditions, logic);
  }, [conditions, logic, onApply]);

  const handleClear = useCallback(() => {
    setConditions([
      { id: genId(), metric: "ochiai", operator: ">=", value: 50 },
    ]);
    setLogic("AND");
    onClear();
  }, [onClear]);

  return (
    <div className="sidebar-section">
      <div
        onClick={() => setExpanded(!expanded)}
        style={{
          display: "flex",
          alignItems: "center",
          justifyContent: "space-between",
          cursor: "pointer",
          marginBottom: expanded ? 12 : 0,
        }}
      >
        <label style={{ cursor: "pointer", marginBottom: 0 }}>Advanced Filter</label>
        <span style={{ fontSize: 12, color: "var(--text-secondary)" }}>
          {expanded ? "[-]" : "[+]"}
        </span>
      </div>

      {expanded && (
        <div style={{ display: "flex", flexDirection: "column", gap: 12 }}>
          {/* Presets */}
          <div>
            <div style={{ fontSize: 11, color: "var(--text-secondary)", marginBottom: 4 }}>
              Quick Presets
            </div>
            <div style={{ display: "flex", gap: 4, flexWrap: "wrap" }}>
              {presets.map((preset) => (
                <button
                  key={preset.name}
                  onClick={() => applyPreset(preset)}
                  style={{
                    padding: "4px 8px",
                    fontSize: 10,
                    border: "1px solid var(--border-color)",
                    borderRadius: 12,
                    background: "white",
                    cursor: "pointer",
                  }}
                >
                  {preset.name}
                </button>
              ))}
            </div>
          </div>

          {/* Logic toggle */}
          <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
            <span style={{ fontSize: 11, color: "var(--text-secondary)" }}>Combine with:</span>
            <div
              style={{
                display: "flex",
                background: "var(--bg-primary)",
                borderRadius: 4,
                border: "1px solid var(--border-color)",
              }}
            >
              <button
                onClick={() => setLogic("AND")}
                style={{
                  padding: "4px 12px",
                  fontSize: 11,
                  border: "none",
                  background: logic === "AND" ? "var(--accent)" : "transparent",
                  color: logic === "AND" ? "white" : "var(--text-primary)",
                  borderRadius: "3px 0 0 3px",
                  cursor: "pointer",
                }}
              >
                AND
              </button>
              <button
                onClick={() => setLogic("OR")}
                style={{
                  padding: "4px 12px",
                  fontSize: 11,
                  border: "none",
                  background: logic === "OR" ? "var(--accent)" : "transparent",
                  color: logic === "OR" ? "white" : "var(--text-primary)",
                  borderRadius: "0 3px 3px 0",
                  cursor: "pointer",
                }}
              >
                OR
              </button>
            </div>
          </div>

          {/* Filter conditions */}
          <div style={{ display: "flex", flexDirection: "column", gap: 8 }}>
            {conditions.map((condition, index) => (
              <div key={condition.id}>
                {index > 0 && (
                  <div
                    style={{
                      textAlign: "center",
                      fontSize: 10,
                      color: "var(--accent)",
                      fontWeight: 600,
                      margin: "4px 0",
                    }}
                  >
                    {logic}
                  </div>
                )}
                <FilterConditionRow
                  condition={condition}
                  onChange={(updated) => updateCondition(condition.id, updated)}
                  onRemove={() => removeCondition(condition.id)}
                  canRemove={conditions.length > 1}
                />
              </div>
            ))}
          </div>

          {/* Add condition button */}
          <button
            onClick={addCondition}
            style={{
              padding: "6px 12px",
              fontSize: 11,
              border: "1px dashed var(--border-color)",
              borderRadius: 4,
              background: "transparent",
              color: "var(--text-secondary)",
              cursor: "pointer",
            }}
          >
            + Add Condition
          </button>

          {/* Action buttons */}
          <div style={{ display: "flex", gap: 8 }}>
            <button
              onClick={handleApply}
              disabled={!state.info}
              style={{
                flex: 1,
                padding: "8px 16px",
                background: "var(--accent)",
                color: "white",
                border: "none",
                borderRadius: 4,
                fontSize: 12,
                fontWeight: 500,
                cursor: "pointer",
              }}
            >
              Apply Filter
            </button>
            <button
              onClick={handleClear}
              style={{
                padding: "8px 16px",
                background: "transparent",
                color: "var(--text-secondary)",
                border: "1px solid var(--border-color)",
                borderRadius: 4,
                fontSize: 12,
                cursor: "pointer",
              }}
            >
              Clear
            </button>
          </div>
        </div>
      )}
    </div>
  );
}
