import { useState, useRef, useEffect } from "react";
import { useDashboard } from "../../state/context";

export default function CutoffSlider() {
  const { state, dispatch } = useDashboard();
  const [local, setLocal] = useState(state.cutoff);
  const debounceRef = useRef<ReturnType<typeof setTimeout>>();

  useEffect(() => setLocal(state.cutoff), [state.cutoff]);

  function handleChange(val: number) {
    setLocal(val);
    clearTimeout(debounceRef.current);
    debounceRef.current = setTimeout(() => {
      dispatch({ type: "SET_CUTOFF", cutoff: val });
    }, 400);
  }

  // Get histogram for current metric if available
  const histogram = state.info?.statistics?.[state.metric]?.histogram;
  const maxBucket = histogram ? Math.max(...histogram) : 1;

  return (
    <div className="sidebar-section">
      <label>Cutoff ({state.metric})</label>
      {histogram && (
        <div className="cutoff-histogram" style={{ display: "flex", alignItems: "flex-end", height: "24px", gap: "1px", marginBottom: "4px" }}>
          {histogram.slice(0, 20).map((count, i) => (
            <div
              key={i}
              style={{
                flex: 1,
                height: `${(count / maxBucket) * 100}%`,
                background: i * 5 < local ? "var(--accent)" : "var(--border-color)",
                minHeight: count > 0 ? "2px" : 0,
              }}
              title={`${i * 5}-${(i + 1) * 5}%: ${count.toLocaleString()} pairs`}
            />
          ))}
        </div>
      )}
      <div className="cutoff-slider">
        <input
          type="range"
          min={0}
          max={100}
          step={1}
          value={local}
          onChange={(e) => handleChange(Number(e.target.value))}
          aria-label="Similarity cutoff threshold"
          aria-valuemin={0}
          aria-valuemax={100}
          aria-valuenow={local}
        />
        <span className="cutoff-value">{local}</span>
      </div>
    </div>
  );
}
