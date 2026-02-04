import { useDashboard } from "../../state/context";
import { LayoutAlgorithm } from "../../state/types";

const LAYOUT_OPTIONS: { value: LayoutAlgorithm; label: string; description: string }[] = [
  { value: "fr", label: "Fruchterman-Reingold", description: "Force-directed layout (default)" },
  { value: "drl", label: "DrL", description: "Distributed Recursive Layout for large graphs" },
  { value: "kk", label: "Kamada-Kawai", description: "Spring-based energy minimization" },
  { value: "circle", label: "Circle", description: "Nodes arranged in a circle" },
  { value: "grid", label: "Grid", description: "Nodes arranged in a grid" },
];

export default function Toolbar() {
  const { state, dispatch } = useDashboard();
  const views = ["network", "table", "stats"] as const;

  const handleLayoutChange = (algorithm: LayoutAlgorithm) => {
    dispatch({ type: "SET_LAYOUT_ALGORITHM", algorithm });
  };

  return (
    <div className="toolbar">
      <div className="view-tabs">
        {views.map((v) => (
          <button
            key={v}
            className={`view-tab ${state.activeView === v ? "active" : ""}`}
            onClick={() => dispatch({ type: "SET_ACTIVE_VIEW", view: v })}
          >
            {v.charAt(0).toUpperCase() + v.slice(1)}
          </button>
        ))}
      </div>

      {/* Layout selector - only shown for network view */}
      {state.activeView === "network" && (
        <div style={{ display: "flex", alignItems: "center", gap: 8, marginLeft: 16 }}>
          <label style={{ fontSize: 12, color: "var(--text-secondary)" }}>Layout:</label>
          <select
            value={state.layoutAlgorithm}
            onChange={(e) => handleLayoutChange(e.target.value as LayoutAlgorithm)}
            style={{
              padding: "4px 8px",
              fontSize: 12,
              border: "1px solid var(--border-color)",
              borderRadius: 4,
              background: "white",
              cursor: "pointer",
            }}
            title={LAYOUT_OPTIONS.find((o) => o.value === state.layoutAlgorithm)?.description}
          >
            {LAYOUT_OPTIONS.map((opt) => (
              <option key={opt.value} value={opt.value} title={opt.description}>
                {opt.label}
              </option>
            ))}
          </select>
        </div>
      )}

      <div className="spacer" />
      {state.graphData && (
        <span className="graph-info">
          {state.graphData.meta.returned_nodes.toLocaleString()} nodes &middot;{" "}
          {state.graphData.meta.returned_edges.toLocaleString()} edges
          {state.focusGroup && (
            <> &middot; focused on &ldquo;{state.focusGroup}&rdquo;</>
          )}
        </span>
      )}
      <button
        className="view-tab"
        style={{ marginLeft: 8 }}
        onClick={() => dispatch({ type: "TOGGLE_QUERY_PANEL" })}
      >
        Query
      </button>
    </div>
  );
}
