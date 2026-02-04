import { useDashboard } from "../../state/context";
import NodeDetail from "./NodeDetail";
import EdgeDetail from "./EdgeDetail";

export default function DetailPanel() {
  const { state, dispatch } = useDashboard();

  function handleClose() {
    dispatch({ type: "SELECT_NODE", node: null });
    dispatch({ type: "SELECT_EDGE", edge: null });
  }

  return (
    <div className="detail-panel">
      <div className="detail-panel-header">
        <h3>{state.selectedNode ? "Node Details" : "Edge Details"}</h3>
        <button className="detail-close" onClick={handleClose}>
          &times;
        </button>
      </div>
      {state.selectedNode && <NodeDetail />}
      {state.selectedEdge && <EdgeDetail />}
    </div>
  );
}
