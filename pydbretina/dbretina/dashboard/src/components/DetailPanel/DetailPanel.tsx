import { useDashboard } from "../../state/context";
import NodeDetail from "./NodeDetail";
import EdgeDetail from "./EdgeDetail";
import ComparePanel from "./ComparePanel";

export default function DetailPanel() {
  const { state, dispatch } = useDashboard();

  const isComparing = state.compareState.active;
  const hasContent = isComparing || !!(state.selectedNode || state.selectedEdge);
  const isOpen = state.detailPanelOpen;

  function handleClose() {
    if (isComparing) {
      dispatch({ type: "EXIT_COMPARE_MODE" });
      dispatch({ type: "SET_PATH_RESULT", result: null });
    } else {
      dispatch({ type: "SELECT_NODE", node: null });
      dispatch({ type: "SELECT_EDGE", edge: null });
    }
    dispatch({ type: "SET_DETAIL_PANEL_OPEN", open: false });
  }

  function handleToggle() {
    dispatch({ type: "TOGGLE_DETAIL_PANEL" });
  }

  return (
    <div className={`detail-panel-wrapper ${isOpen ? "open" : "collapsed"}`}>
      {/* Toggle strip — always visible */}
      <button
        className="detail-panel-toggle"
        onClick={handleToggle}
        title={isOpen ? "Collapse panel" : "Expand panel"}
      >
        <svg width="8" height="14" viewBox="0 0 8 14" fill="none">
          {isOpen ? (
            <path d="M1 1l6 6-6 6" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
          ) : (
            <path d="M7 1L1 7l6 6" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
          )}
        </svg>
        {!isOpen && hasContent && <span className="detail-panel-badge" />}
      </button>

      {/* Panel content */}
      <div className="detail-panel">
        {isComparing ? (
          <ComparePanel />
        ) : (
          <>
            <div className="detail-panel-header">
              <h3>
                {state.selectedNode
                  ? "Node Details"
                  : state.selectedEdge
                    ? "Edge Details"
                    : "Details"}
              </h3>
              <button className="detail-close" onClick={handleClose}>
                &times;
              </button>
            </div>
            {state.selectedNode && <NodeDetail />}
            {state.selectedEdge && <EdgeDetail />}
            {!state.selectedNode && !state.selectedEdge && (
              <div
                style={{
                  padding: "32px 16px",
                  textAlign: "center",
                  color: "var(--text-secondary)",
                  fontSize: 13,
                }}
              >
                <svg
                  width="32"
                  height="32"
                  viewBox="0 0 24 24"
                  fill="none"
                  stroke="currentColor"
                  strokeWidth="1.5"
                  style={{ opacity: 0.4, marginBottom: 12 }}
                >
                  <circle cx="12" cy="12" r="10" />
                  <path d="M12 16v-4M12 8h.01" />
                </svg>
                <p>Click a node or edge in the graph to see details here.</p>
              </div>
            )}
          </>
        )}
      </div>
    </div>
  );
}
