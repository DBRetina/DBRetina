import { useDashboard } from "../../state/context";

export default function DatasetInfo() {
  const { state } = useDashboard();
  const info = state.info;
  if (!info) return null;

  return (
    <div className="sidebar-section">
      <label>Dataset</label>
      <div className="info-grid">
        <div className="info-item">
          <div className="info-value">{info.num_groups.toLocaleString()}</div>
          <div className="info-label">Groups</div>
        </div>
        <div className="info-item">
          <div className="info-value">{info.num_pairs.toLocaleString()}</div>
          <div className="info-label">Pairs</div>
        </div>
      </div>
      {state.graphData && (
        <div className="info-grid" style={{ marginTop: 6 }}>
          <div className="info-item">
            <div className="info-value">
              {state.graphData.meta.returned_nodes.toLocaleString()}
            </div>
            <div className="info-label">Visible nodes</div>
          </div>
          <div className="info-item">
            <div className="info-value">
              {state.graphData.meta.returned_edges.toLocaleString()}
            </div>
            <div className="info-label">Visible edges</div>
          </div>
        </div>
      )}
      {state.graphData?.meta.truncated && (
        <div style={{ fontSize: 11, color: "var(--warning)", marginTop: 4 }}>
          Graph truncated to top-degree nodes
        </div>
      )}
    </div>
  );
}
