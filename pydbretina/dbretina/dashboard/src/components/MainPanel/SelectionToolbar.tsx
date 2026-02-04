import { useDashboard } from "../../state/context";
import { SelectionMode } from "../../state/types";

const SELECTION_MODES: { mode: SelectionMode; label: string; icon: string; title: string }[] = [
  { mode: "click", label: "Click", icon: "👆", title: "Click to select individual nodes" },
  { mode: "lasso", label: "Lasso", icon: "🔲", title: "Draw to select multiple nodes" },
  { mode: "box", label: "Box", icon: "⬜", title: "Drag a rectangle to select nodes" },
];

export default function SelectionToolbar() {
  const { state, dispatch } = useDashboard();
  const { selection, graphData } = state;

  const handleModeChange = (mode: SelectionMode) => {
    dispatch({ type: "SET_SELECTION_MODE", mode });
  };

  const handleClearSelection = () => {
    dispatch({ type: "CLEAR_SELECTION" });
  };

  const handleInvertSelection = () => {
    dispatch({ type: "INVERT_SELECTION" });
  };

  const handleSelectAll = () => {
    if (!graphData) return;
    dispatch({ type: "SET_SELECTED_NODES", nodes: new Set(graphData.nodes.map((n) => n.id)) });
  };

  const selectedCount = selection.selectedNodes.size;

  return (
    <div
      style={{
        display: "flex",
        alignItems: "center",
        gap: 12,
        padding: "8px 12px",
        background: "var(--bg-secondary)",
        borderBottom: "1px solid var(--border-color)",
        fontSize: 12,
      }}
    >
      {/* Selection mode buttons */}
      <div style={{ display: "flex", gap: 4 }}>
        {SELECTION_MODES.map(({ mode, label, icon, title }) => (
          <button
            key={mode}
            onClick={() => handleModeChange(mode)}
            title={title}
            style={{
              padding: "4px 8px",
              border: `1px solid ${selection.mode === mode ? "var(--accent)" : "var(--border-color)"}`,
              borderRadius: 4,
              background: selection.mode === mode ? "var(--accent)" : "white",
              color: selection.mode === mode ? "white" : "var(--text-primary)",
              cursor: "pointer",
              display: "flex",
              alignItems: "center",
              gap: 4,
              fontSize: 11,
            }}
          >
            <span>{icon}</span>
            <span>{label}</span>
          </button>
        ))}
      </div>

      <div style={{ width: 1, height: 20, background: "var(--border-color)" }} />

      {/* Selection info and actions */}
      <span style={{ color: "var(--text-secondary)" }}>
        {selectedCount > 0 ? `${selectedCount} selected` : "No selection"}
      </span>

      {selectedCount > 0 && (
        <>
          <button
            onClick={handleClearSelection}
            style={{
              padding: "4px 8px",
              border: "1px solid var(--border-color)",
              borderRadius: 4,
              background: "white",
              cursor: "pointer",
              fontSize: 11,
            }}
          >
            Clear
          </button>
          <button
            onClick={handleInvertSelection}
            style={{
              padding: "4px 8px",
              border: "1px solid var(--border-color)",
              borderRadius: 4,
              background: "white",
              cursor: "pointer",
              fontSize: 11,
            }}
          >
            Invert
          </button>
        </>
      )}

      {selectedCount === 0 && graphData && (
        <button
          onClick={handleSelectAll}
          style={{
            padding: "4px 8px",
            border: "1px solid var(--border-color)",
            borderRadius: 4,
            background: "white",
            cursor: "pointer",
            fontSize: 11,
          }}
        >
          Select All
        </button>
      )}

      <div style={{ flex: 1 }} />

      {/* Selection actions */}
      {selectedCount > 0 && (
        <div style={{ display: "flex", gap: 4 }}>
          <button
            onClick={() => {
              // Highlight selected nodes
              dispatch({ type: "SET_HIGHLIGHT_NODES", nodes: selection.selectedNodes });
            }}
            style={{
              padding: "4px 8px",
              border: "1px solid var(--accent)",
              borderRadius: 4,
              background: "white",
              color: "var(--accent)",
              cursor: "pointer",
              fontSize: 11,
            }}
          >
            Highlight
          </button>
          <button
            onClick={() => {
              // Export selected node IDs
              const nodeIds = Array.from(selection.selectedNodes);
              const blob = new Blob([JSON.stringify(nodeIds, null, 2)], { type: "application/json" });
              const url = URL.createObjectURL(blob);
              const a = document.createElement("a");
              a.href = url;
              a.download = "selected_nodes.json";
              a.click();
              URL.revokeObjectURL(url);
            }}
            style={{
              padding: "4px 8px",
              border: "1px solid var(--border-color)",
              borderRadius: 4,
              background: "white",
              cursor: "pointer",
              fontSize: 11,
            }}
          >
            Export
          </button>
        </div>
      )}
    </div>
  );
}
