import { useMemo } from "react";
import { useDashboard } from "../../state/context";
import { LayoutAlgorithm } from "../../state/types";

const LAYOUT_OPTIONS: { value: LayoutAlgorithm; label: string; description: string }[] = [
  { value: "force", label: "Force-Directed", description: "Live force simulation (interactive)" },
  { value: "fr", label: "Fruchterman-Reingold", description: "Classic FR layout (server-computed)" },
  { value: "fa2", label: "ForceAtlas2", description: "ForceAtlas2 layout (good for communities)" },
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

  // Visible node/edge counts, accounting for the client-side node/edge filters
  // (Advanced Filter). When no filter is active these equal the server-returned
  // counts. Mirrors the compose logic in useForceGraphData.
  const visible = useMemo(() => {
    const gd = state.graphData;
    if (!gd) return { nodes: 0, edges: 0, filtered: false };
    const { nodeFilter, edgeFilter } = state;
    if (!nodeFilter && !edgeFilter) {
      return {
        nodes: gd.meta.returned_nodes,
        edges: gd.meta.returned_edges,
        filtered: false,
      };
    }
    let edges = gd.edges;
    if (edgeFilter) edges = edges.filter((e) => edgeFilter.has(`${e.source}|${e.target}`));
    if (nodeFilter)
      edges = edges.filter((e) => nodeFilter.has(e.source) && nodeFilter.has(e.target));
    let nodeIds: Set<string>;
    if (edgeFilter) {
      nodeIds = new Set<string>();
      for (const e of edges) {
        nodeIds.add(e.source);
        nodeIds.add(e.target);
      }
      if (nodeFilter) nodeIds = new Set([...nodeIds].filter((id) => nodeFilter.has(id)));
    } else {
      // node-filter only: count filtered nodes directly
      nodeIds = new Set(gd.nodes.filter((n) => nodeFilter!.has(n.id)).map((n) => n.id));
    }
    return { nodes: nodeIds.size, edges: edges.length, filtered: true };
  }, [state.graphData, state.nodeFilter, state.edgeFilter]);

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
          <span data-testid="visible-nodes">{visible.nodes.toLocaleString()}</span> nodes &middot;{" "}
          <span data-testid="visible-edges">{visible.edges.toLocaleString()}</span> edges
          {visible.filtered && (
            <span data-testid="filter-active"> &middot; filtered</span>
          )}
          {state.focusGroup && (
            <> &middot; focused on &ldquo;{state.focusGroup}&rdquo;</>
          )}
        </span>
      )}
      <button
        className={`view-tab ${state.queryPanelOpen ? "active" : ""}`}
        style={{ marginLeft: 8 }}
        onClick={() => dispatch({ type: "TOGGLE_QUERY_PANEL" })}
      >
        Query
      </button>
    </div>
  );
}
