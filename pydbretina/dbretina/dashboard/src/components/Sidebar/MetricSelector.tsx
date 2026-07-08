import { useDashboard } from "../../state/context";

export default function MetricSelector() {
  const { state, dispatch } = useDashboard();
  const metrics = state.info?.valid_metrics ?? [];

  return (
    <div className="sidebar-section">
      <label>Metric</label>
      <select
        value={state.metric}
        onChange={(e) => dispatch({ type: "SET_METRIC", metric: e.target.value })}
        aria-label="Similarity metric"
      >
        {metrics.map((m) => (
          <option key={m} value={m}>
            {m}
          </option>
        ))}
      </select>
    </div>
  );
}
