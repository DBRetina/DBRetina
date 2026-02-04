import { useState, useEffect } from "react";
import { useDashboard } from "../../state/context";

export default function NodeDetail() {
  const { state, dispatch } = useDashboard();
  const node = state.selectedNode;
  if (!node) return null;

  // Find connected edges
  const edges = state.graphData?.edges.filter(
    (e) => e.source === node.id || e.target === node.id
  ) ?? [];

  // Top neighbors by weight
  const neighbors = edges
    .map((e) => {
      const otherId = e.source === node.id ? e.target : e.source;
      const otherNode = state.graphData?.nodes.find((n) => n.id === otherId);
      return { name: otherNode?.label ?? otherId, weight: e.weight };
    })
    .sort((a, b) => b.weight - a.weight)
    .slice(0, 10);

  function handleFocus() {
    dispatch({ type: "SET_FOCUS_GROUP", group: node!.label });
    dispatch({ type: "SELECT_NODE", node: null });
  }

  return (
    <>
      <div className="detail-section">
        <h4>Properties</h4>
        <div className="detail-row">
          <span className="detail-label">Name</span>
          <span className="detail-value" style={{ fontSize: 11, maxWidth: 180, overflow: "hidden", textOverflow: "ellipsis" }}>
            {node.label}
          </span>
        </div>
        <div className="detail-row">
          <span className="detail-label">Degree</span>
          <span className="detail-value">{node.degree}</span>
        </div>
        <div className="detail-row">
          <span className="detail-label">Community</span>
          <span className="detail-value">{node.community}</span>
        </div>
        <div className="detail-row">
          <span className="detail-label">PageRank</span>
          <span className="detail-value">{node.pagerank.toFixed(6)}</span>
        </div>
        <button
          onClick={handleFocus}
          style={{
            marginTop: 8,
            padding: "5px 12px",
            fontSize: 11,
            background: "var(--accent)",
            color: "white",
            border: "none",
            borderRadius: 4,
            cursor: "pointer",
          }}
        >
          Focus neighborhood
        </button>
      </div>
      {neighbors.length > 0 && (
        <div className="detail-section">
          <h4>Top Neighbors ({state.metric})</h4>
          {neighbors.map((n, i) => (
            <div key={i} className="detail-row">
              <span className="detail-label" style={{ maxWidth: 180, overflow: "hidden", textOverflow: "ellipsis", whiteSpace: "nowrap" }}>
                {n.name}
              </span>
              <span className="detail-value">{n.weight.toFixed(1)}</span>
            </div>
          ))}
        </div>
      )}
    </>
  );
}
