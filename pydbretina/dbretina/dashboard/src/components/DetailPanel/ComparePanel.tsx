import { useState, useCallback } from "react";
import { useDashboard } from "../../state/context";
import {
  fetchSharedFeatures,
  explainPair,
  fetchShortestPath,
  fetchMetricProfile,
  MetricProfileResult,
} from "../../api";
import CollapsibleSection from "./CollapsibleSection";

type ExplainMethod = "hypergraph" | "edge_weighted" | "projection";

export default function ComparePanel() {
  const { state, dispatch } = useDashboard();
  // nodeA is always set when ComparePanel is rendered (compareState.active === true)
  const nodeA = state.compareState.nodeA!;
  const nodeB = state.compareState.nodeB;

  // Shared genes state
  const [sharedGenes, setSharedGenes] = useState<string[] | null>(null);
  const [sharedCount, setSharedCount] = useState(0);
  const [sharedLoading, setSharedLoading] = useState(false);
  const [sharedError, setSharedError] = useState<string | null>(null);

  // Gene importance state
  const [explainMethod, setExplainMethod] = useState<ExplainMethod>("hypergraph");
  const [importantGenes, setImportantGenes] = useState<
    Array<{ gene: string; score: number; specificity?: number }> | null
  >(null);
  const [importLoading, setImportLoading] = useState(false);
  const [importError, setImportError] = useState<string | null>(null);

  // Shortest path state
  const [pathNodes, setPathNodes] = useState<string[] | null>(null);
  const [pathLoading, setPathLoading] = useState(false);
  const [pathError, setPathError] = useState<string | null>(null);

  // Metric comparison state
  const [profileA, setProfileA] = useState<MetricProfileResult | null>(null);
  const [profileB, setProfileB] = useState<MetricProfileResult | null>(null);
  const [metricLoading, setMetricLoading] = useState(false);
  const [metricError, setMetricError] = useState<string | null>(null);

  // Auto-load shared genes when both nodes are set
  const fetchShared = useCallback(async () => {
    if (!nodeA || !nodeB) return;
    setSharedLoading(true);
    setSharedError(null);
    try {
      const res = await fetchSharedFeatures(nodeA.label, nodeB.label);
      setSharedGenes(res.features.slice(0, 50));
      setSharedCount(res.count);
    } catch (err) {
      setSharedError(err instanceof Error ? err.message : "Failed to load");
      setSharedGenes(null);
    } finally {
      setSharedLoading(false);
    }
  }, [nodeA, nodeB]);

  // Gene importance
  const fetchImportance = useCallback(
    async (method?: ExplainMethod) => {
      if (!nodeA || !nodeB) return;
      const m = method ?? explainMethod;
      setImportLoading(true);
      setImportError(null);
      try {
        const res = await explainPair(nodeA.label, nodeB.label, m);
        setImportantGenes(res.genes);
      } catch (err) {
        setImportError(err instanceof Error ? err.message : "Failed to load");
        setImportantGenes(null);
      } finally {
        setImportLoading(false);
      }
    },
    [nodeA, nodeB, explainMethod]
  );

  // Shortest path
  const fetchPath = useCallback(async () => {
    if (!nodeA || !nodeB) return;
    setPathLoading(true);
    setPathError(null);
    try {
      const res = await fetchShortestPath(
        nodeA.label,
        nodeB.label,
        state.metric,
        state.cutoff
      );
      setPathNodes(res.connected ? res.path_nodes : []);
      // Highlight path nodes in graph
      if (res.connected && res.path_nodes.length > 0) {
        dispatch({ type: "SET_PATH_RESULT", result: res });
      }
    } catch (err) {
      setPathError(err instanceof Error ? err.message : "Failed to load");
      setPathNodes(null);
    } finally {
      setPathLoading(false);
    }
  }, [nodeA, nodeB, state.metric, state.cutoff, dispatch]);

  // Metric comparison
  const fetchMetrics = useCallback(async () => {
    if (!nodeA || !nodeB) return;
    setMetricLoading(true);
    setMetricError(null);
    try {
      const [a, b] = await Promise.all([
        fetchMetricProfile(nodeA.label),
        fetchMetricProfile(nodeB.label),
      ]);
      setProfileA(a);
      setProfileB(b);
    } catch (err) {
      setMetricError(err instanceof Error ? err.message : "Failed to load");
      setProfileA(null);
      setProfileB(null);
    } finally {
      setMetricLoading(false);
    }
  }, [nodeA, nodeB]);

  function handleClose() {
    dispatch({ type: "EXIT_COMPARE_MODE" });
    dispatch({ type: "SET_PATH_RESULT", result: null });
  }

  // Waiting for node B
  if (!nodeB) {
    return (
      <>
        <div className="detail-panel-header">
          <h3>Compare Mode</h3>
          <button className="detail-close" onClick={handleClose}>
            &times;
          </button>
        </div>
        <div className="detail-section">
          <div style={{ fontSize: 13, marginBottom: 8 }}>
            Click a second node to compare with:
          </div>
          <div
            style={{
              padding: 12,
              background: "var(--bg-primary)",
              borderRadius: 6,
              border: "1px solid var(--accent)",
              fontSize: 14,
              fontWeight: 600,
              color: "var(--accent)",
            }}
          >
            {nodeA?.label}
          </div>
          <button
            onClick={handleClose}
            style={{
              marginTop: 12,
              width: "100%",
              padding: "8px 16px",
              background: "white",
              color: "var(--text-primary)",
              border: "1px solid var(--border-color)",
              borderRadius: 4,
              fontSize: 12,
              cursor: "pointer",
            }}
          >
            Cancel
          </button>
        </div>
      </>
    );
  }

  // Both nodes selected — show comparison sections
  return (
    <>
      <div className="detail-panel-header">
        <h3>Compare</h3>
        <button className="detail-close" onClick={handleClose}>
          &times;
        </button>
      </div>

      {/* Compared nodes */}
      <div className="detail-section">
        <div
          style={{
            display: "flex",
            alignItems: "center",
            gap: 8,
            fontSize: 12,
          }}
        >
          <div
            style={{
              flex: 1,
              padding: 8,
              background: "var(--bg-primary)",
              borderRadius: 4,
              fontWeight: 600,
              textAlign: "center",
              overflow: "hidden",
              textOverflow: "ellipsis",
              whiteSpace: "nowrap",
            }}
          >
            {nodeA.label}
          </div>
          <span style={{ color: "var(--text-secondary)", fontSize: 16 }}>&harr;</span>
          <div
            style={{
              flex: 1,
              padding: 8,
              background: "var(--bg-primary)",
              borderRadius: 4,
              fontWeight: 600,
              textAlign: "center",
              overflow: "hidden",
              textOverflow: "ellipsis",
              whiteSpace: "nowrap",
            }}
          >
            {nodeB.label}
          </div>
        </div>
      </div>

      {/* Shared Genes */}
      {state.info?.has_genes && (
        <CollapsibleSection
          title="Shared Genes"
          defaultExpanded
          loading={sharedLoading}
          error={sharedError}
          onExpand={fetchShared}
          onRetry={fetchShared}
          badge={sharedCount || undefined}
        >
          {sharedGenes && sharedGenes.length > 0 ? (
            <div className="gene-list">
              {sharedGenes.map((g) => (
                <span key={g} className="gene-tag">
                  {g}
                </span>
              ))}
              {sharedCount > 50 && (
                <span className="gene-tag" style={{ color: "var(--text-secondary)" }}>
                  +{sharedCount - 50} more
                </span>
              )}
            </div>
          ) : sharedGenes ? (
            <div style={{ fontSize: 11, color: "var(--text-secondary)" }}>
              No shared genes found
            </div>
          ) : null}
        </CollapsibleSection>
      )}

      {/* Gene Importance */}
      {state.info?.has_genes && (
        <CollapsibleSection
          title="Gene Importance"
          loading={importLoading}
          error={importError}
          onExpand={() => fetchImportance()}
          onRetry={() => fetchImportance()}
          badge={importantGenes ? importantGenes.length : undefined}
        >
          {/* Method toggles */}
          <div
            style={{
              display: "flex",
              gap: 2,
              marginBottom: 8,
              background: "var(--bg-primary)",
              borderRadius: 4,
              padding: 2,
            }}
          >
            {(
              [
                ["hypergraph", "TF-IDF"],
                ["edge_weighted", "Edge-Wtd"],
                ["projection", "PageRank"],
              ] as [ExplainMethod, string][]
            ).map(([method, label]) => (
              <button
                key={method}
                onClick={() => {
                  setExplainMethod(method);
                  if (importantGenes) fetchImportance(method);
                }}
                style={{
                  flex: 1,
                  padding: "4px 8px",
                  border: "none",
                  borderRadius: 3,
                  fontSize: 10,
                  fontWeight: 500,
                  cursor: "pointer",
                  background: explainMethod === method ? "var(--accent)" : "transparent",
                  color: explainMethod === method ? "white" : "var(--text-secondary)",
                }}
              >
                {label}
              </button>
            ))}
          </div>

          {importantGenes && importantGenes.length > 0 && (
            <div
              style={{
                maxHeight: 200,
                overflowY: "auto",
                border: "1px solid var(--border-color)",
                borderRadius: 4,
              }}
            >
              <table style={{ width: "100%", fontSize: 11, borderCollapse: "collapse" }}>
                <thead>
                  <tr style={{ background: "var(--bg-primary)" }}>
                    <th style={{ padding: "5px 8px", textAlign: "left" }}>Gene</th>
                    <th style={{ padding: "5px 8px", textAlign: "right" }}>Score</th>
                  </tr>
                </thead>
                <tbody>
                  {importantGenes.slice(0, 30).map((g, i) => (
                    <tr
                      key={g.gene}
                      style={{
                        borderBottom:
                          i < Math.min(importantGenes.length, 30) - 1
                            ? "1px solid var(--border-color)"
                            : "none",
                      }}
                    >
                      <td
                        style={{
                          padding: "5px 8px",
                          fontWeight: i < 3 ? 600 : 400,
                          color: i < 3 ? "var(--accent)" : "inherit",
                        }}
                      >
                        {g.gene}
                        {g.specificity !== undefined && g.specificity > 0.8 && (
                          <span
                            style={{
                              marginLeft: 4,
                              fontSize: 9,
                              padding: "1px 4px",
                              background: "#d4edda",
                              borderRadius: 2,
                              color: "#155724",
                            }}
                          >
                            specific
                          </span>
                        )}
                      </td>
                      <td
                        style={{
                          padding: "5px 8px",
                          textAlign: "right",
                          fontFamily: "var(--font-mono)",
                          fontSize: 10,
                        }}
                      >
                        {g.score.toFixed(4)}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}
        </CollapsibleSection>
      )}

      {/* Shortest Path */}
      <CollapsibleSection
        title="Shortest Path"
        loading={pathLoading}
        error={pathError}
        onExpand={fetchPath}
        onRetry={fetchPath}
        badge={pathNodes ? pathNodes.length : undefined}
      >
        {pathNodes !== null && (
          pathNodes.length > 0 ? (
            <div style={{ fontSize: 12 }}>
              <div
                style={{
                  display: "flex",
                  flexWrap: "wrap",
                  alignItems: "center",
                  gap: 4,
                }}
              >
                {pathNodes.map((name, i) => (
                  <span key={i} style={{ display: "flex", alignItems: "center", gap: 4 }}>
                    <span
                      style={{
                        padding: "3px 8px",
                        background:
                          name === nodeA.label || name === nodeB.label
                            ? "var(--accent)"
                            : "var(--bg-primary)",
                        color:
                          name === nodeA.label || name === nodeB.label
                            ? "white"
                            : "var(--text-primary)",
                        borderRadius: 4,
                        fontSize: 11,
                        fontWeight: 500,
                      }}
                    >
                      {name}
                    </span>
                    {i < pathNodes.length - 1 && (
                      <span style={{ color: "var(--text-secondary)", fontSize: 10 }}>
                        &rarr;
                      </span>
                    )}
                  </span>
                ))}
              </div>
              <div style={{ fontSize: 10, color: "var(--text-secondary)", marginTop: 6 }}>
                Path length: {pathNodes.length - 1} hops
              </div>
            </div>
          ) : (
            <div style={{ fontSize: 11, color: "var(--text-secondary)" }}>
              No path found — nodes are not connected at current cutoff
            </div>
          )
        )}
      </CollapsibleSection>

      {/* Metric Comparison */}
      <CollapsibleSection
        title="Metric Comparison"
        loading={metricLoading}
        error={metricError}
        onExpand={fetchMetrics}
        onRetry={fetchMetrics}
      >
        {profileA && profileB && (
          <div style={{ display: "flex", flexDirection: "column", gap: 8 }}>
            {(() => {
              const metricsA = new Map(profileA.metrics.map((m) => [m.metric, m]));
              const metricsB = new Map(profileB.metrics.map((m) => [m.metric, m]));
              const allMetrics = [...new Set([...metricsA.keys(), ...metricsB.keys()])];
              const maxVal = Math.max(
                ...profileA.metrics.map((m) => m.avg),
                ...profileB.metrics.map((m) => m.avg),
                0.001
              );

              return allMetrics.map((metric) => {
                const a = metricsA.get(metric);
                const b = metricsB.get(metric);
                return (
                  <div key={metric} style={{ fontSize: 11 }}>
                    <div style={{ fontWeight: 500, marginBottom: 3 }}>{metric}</div>
                    {/* Node A bar */}
                    <div style={{ display: "flex", alignItems: "center", gap: 6, marginBottom: 2 }}>
                      <span style={{ width: 12, fontSize: 9, color: "var(--text-secondary)" }}>A</span>
                      <div
                        style={{
                          flex: 1,
                          height: 5,
                          background: "var(--bg-primary)",
                          borderRadius: 3,
                          overflow: "hidden",
                        }}
                      >
                        <div
                          style={{
                            height: "100%",
                            width: `${a ? (a.avg / maxVal) * 100 : 0}%`,
                            background: "var(--accent)",
                            borderRadius: 3,
                          }}
                        />
                      </div>
                      <span style={{ fontSize: 9, fontFamily: "var(--font-mono)", minWidth: 36, textAlign: "right" }}>
                        {a?.avg.toFixed(2) ?? "—"}
                      </span>
                    </div>
                    {/* Node B bar */}
                    <div style={{ display: "flex", alignItems: "center", gap: 6 }}>
                      <span style={{ width: 12, fontSize: 9, color: "var(--text-secondary)" }}>B</span>
                      <div
                        style={{
                          flex: 1,
                          height: 5,
                          background: "var(--bg-primary)",
                          borderRadius: 3,
                          overflow: "hidden",
                        }}
                      >
                        <div
                          style={{
                            height: "100%",
                            width: `${b ? (b.avg / maxVal) * 100 : 0}%`,
                            background: "var(--success)",
                            borderRadius: 3,
                          }}
                        />
                      </div>
                      <span style={{ fontSize: 9, fontFamily: "var(--font-mono)", minWidth: 36, textAlign: "right" }}>
                        {b?.avg.toFixed(2) ?? "—"}
                      </span>
                    </div>
                  </div>
                );
              });
            })()}
            <div style={{ display: "flex", gap: 16, fontSize: 10, color: "var(--text-secondary)", marginTop: 4 }}>
              <span>
                <span style={{ display: "inline-block", width: 8, height: 8, background: "var(--accent)", borderRadius: 2, marginRight: 4 }} />
                {nodeA.label}
              </span>
              <span>
                <span style={{ display: "inline-block", width: 8, height: 8, background: "var(--success)", borderRadius: 2, marginRight: 4 }} />
                {nodeB.label}
              </span>
            </div>
          </div>
        )}
      </CollapsibleSection>
    </>
  );
}
