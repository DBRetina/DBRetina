import { useDashboard } from "../../state/context";
import ForceGraphRenderer from "./ForceGraphRenderer";
import SelectionToolbar from "./SelectionToolbar";
import { EmptyState } from "../Common";

export default function NetworkView() {
  const { state } = useDashboard();

  // Check if there's an error state
  const graphError = state.errors.graph;
  if (graphError && !state.graphData) {
    return (
      <div className="graph-container">
        <EmptyState
          type="error"
          title="Failed to load graph"
          message={graphError.detail}
          action={
            graphError.context?.suggestion
              ? { label: "Learn more", onClick: () => {} }
              : undefined
          }
        />
      </div>
    );
  }

  // Check for loading state - show nothing, loading overlay handled by App
  if (state.loadingProgress.graph?.active) {
    return <div className="graph-container" />;
  }

  // Check for empty data
  if (!state.graphData || state.graphData.nodes.length === 0) {
    return (
      <div className="graph-container">
        <EmptyState
          type="no-data"
          title="No graph data"
          message={
            state.focusGroup
              ? `No connections found for "${state.focusGroup}" at this cutoff.`
              : "No pairs meet the current filter criteria. Try lowering the cutoff value."
          }
        />
      </div>
    );
  }

  return (
    <div style={{ display: "flex", flexDirection: "column", height: "100%", width: "100%" }}>
      <SelectionToolbar />
      <div style={{ flex: 1, position: "relative", minHeight: 0 }}>
        <ForceGraphRenderer />
      </div>
    </div>
  );
}
