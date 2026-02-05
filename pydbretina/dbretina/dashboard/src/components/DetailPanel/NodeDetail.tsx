import { useState, useEffect, useCallback, useRef } from "react";
import { useDashboard } from "../../state/context";
import {
  fetchGroupGenes,
  fetchHubGenes,
  fetchMetricProfile,
  analyzeClusterGenes,
  MetricProfileResult,
} from "../../api";
import CollapsibleSection from "./CollapsibleSection";

type HubMethod = "hypergraph" | "edge_weighted" | "projection";

export default function NodeDetail() {
  const { state, dispatch } = useDashboard();
  const node = state.selectedNode;

  // Gene state
  const [genes, setGenes] = useState<string[]>([]);
  const [geneCount, setGeneCount] = useState<number>(0);
  const [loadingGenes, setLoadingGenes] = useState(false);

  // Hub genes state
  const [hubMethod, setHubMethod] = useState<HubMethod>("hypergraph");
  const [hubGenes, setHubGenes] = useState<
    Array<{ gene: string; score: number; specificity?: number }> | null
  >(null);
  const [hubLoading, setHubLoading] = useState(false);
  const [hubError, setHubError] = useState<string | null>(null);
  const hubFetchedFor = useRef<string | null>(null);

  // Metric profile state
  const [metricProfile, setMetricProfile] = useState<MetricProfileResult | null>(null);
  const [profileLoading, setProfileLoading] = useState(false);
  const [profileError, setProfileError] = useState<string | null>(null);

  // Neighborhood enrichment state
  const [enrichGenes, setEnrichGenes] = useState<
    Array<{ gene: string; score: number; in_cluster_count: number; cluster_size: number }> | null
  >(null);
  const [enrichLoading, setEnrichLoading] = useState(false);
  const [enrichError, setEnrichError] = useState<string | null>(null);

  // Reset all section states when node changes
  useEffect(() => {
    setHubGenes(null);
    setHubError(null);
    hubFetchedFor.current = null;
    setMetricProfile(null);
    setProfileError(null);
    setEnrichGenes(null);
    setEnrichError(null);
  }, [node?.label]);

  // Fetch genes when node changes
  useEffect(() => {
    if (!node || !state.info?.has_genes) {
      setGenes([]);
      setGeneCount(0);
      return;
    }

    setLoadingGenes(true);
    fetchGroupGenes(node.label, 50)
      .then((res) => {
        setGenes(res.genes);
        setGeneCount(res.count);
      })
      .catch(() => {
        setGenes([]);
        setGeneCount(0);
      })
      .finally(() => setLoadingGenes(false));
  }, [node?.label, state.info?.has_genes]);

  // Precompute neighbors
  const neighbors = (() => {
    if (!node || !state.graphData) return [];
    const idToLabel = new Map<string, string>();
    for (const n of state.graphData.nodes) {
      idToLabel.set(n.id, n.label);
    }
    const result: { name: string; weight: number }[] = [];
    for (const e of state.graphData.edges) {
      let otherId: string | null = null;
      if (e.source === node.id) otherId = e.target;
      else if (e.target === node.id) otherId = e.source;
      if (otherId !== null) {
        result.push({ name: idToLabel.get(otherId) ?? otherId, weight: e.weight });
      }
    }
    result.sort((a, b) => b.weight - a.weight);
    return result.slice(0, 10);
  })();

  // Hub genes fetch
  const fetchHub = useCallback(
    async (method: HubMethod) => {
      if (!node) return;
      setHubLoading(true);
      setHubError(null);
      try {
        const res = await fetchHubGenes(node.label, method, 2, 20);
        setHubGenes(res.genes);
        hubFetchedFor.current = node.label;
      } catch (err) {
        setHubError(err instanceof Error ? err.message : "Failed to load hub genes");
        setHubGenes(null);
      } finally {
        setHubLoading(false);
      }
    },
    [node]
  );

  // Metric profile fetch
  const fetchProfile = useCallback(async () => {
    if (!node) return;
    setProfileLoading(true);
    setProfileError(null);
    try {
      const res = await fetchMetricProfile(node.label);
      setMetricProfile(res);
    } catch (err) {
      setProfileError(err instanceof Error ? err.message : "Failed to load profile");
      setMetricProfile(null);
    } finally {
      setProfileLoading(false);
    }
  }, [node]);

  // Neighborhood enrichment fetch
  const fetchEnrichment = useCallback(async () => {
    if (!node || !state.graphData) return;
    setEnrichLoading(true);
    setEnrichError(null);
    try {
      // Collect node + direct neighbors
      const idToLabel = new Map<string, string>();
      for (const n of state.graphData.nodes) {
        idToLabel.set(n.id, n.label);
      }
      const neighborLabels = new Set<string>();
      neighborLabels.add(node.label);
      for (const e of state.graphData.edges) {
        if (e.source === node.id) {
          const lbl = idToLabel.get(e.target);
          if (lbl) neighborLabels.add(lbl);
        } else if (e.target === node.id) {
          const lbl = idToLabel.get(e.source);
          if (lbl) neighborLabels.add(lbl);
        }
      }
      const res = await analyzeClusterGenes(Array.from(neighborLabels), "hypergraph", 30);
      setEnrichGenes(res.genes);
    } catch (err) {
      setEnrichError(err instanceof Error ? err.message : "Failed to load enrichment");
      setEnrichGenes(null);
    } finally {
      setEnrichLoading(false);
    }
  }, [node, state.graphData]);

  // Re-fetch hub genes when method changes (only if already fetched)
  useEffect(() => {
    if (node && hubFetchedFor.current === node.label) {
      fetchHub(hubMethod);
    }
  }, [hubMethod]);

  if (!node) return null;

  function handleFocus() {
    dispatch({ type: "SET_FOCUS_GROUP", group: node!.label });
    dispatch({ type: "SELECT_NODE", node: null });
  }

  function handleCompare() {
    dispatch({ type: "ENTER_COMPARE_MODE", nodeA: node! });
  }

  return (
    <>
      {/* Properties */}
      <div className="detail-section">
        <h4>Properties</h4>
        <div className="detail-row">
          <span className="detail-label">Name</span>
          <span
            className="detail-value"
            style={{
              fontSize: 11,
              maxWidth: 180,
              overflow: "hidden",
              textOverflow: "ellipsis",
            }}
          >
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
        <div style={{ display: "flex", gap: 6, marginTop: 8 }}>
          <button
            onClick={handleFocus}
            style={{
              flex: 1,
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
          <button
            onClick={handleCompare}
            style={{
              flex: 1,
              padding: "5px 12px",
              fontSize: 11,
              background: "white",
              color: "var(--accent)",
              border: "1px solid var(--accent)",
              borderRadius: 4,
              cursor: "pointer",
            }}
          >
            Compare with...
          </button>
        </div>
      </div>

      {/* Genes section */}
      {state.info?.has_genes && (
        <div className="detail-section">
          <h4>Genes {geneCount > 0 && `(${geneCount})`}</h4>
          {loadingGenes ? (
            <div style={{ fontSize: 12, color: "var(--text-secondary)" }}>
              Loading...
            </div>
          ) : genes.length > 0 ? (
            <div className="gene-list">
              {genes.map((g) => (
                <span key={g} className="gene-tag">
                  {g}
                </span>
              ))}
              {geneCount > 50 && (
                <span
                  className="gene-tag"
                  style={{ color: "var(--text-secondary)" }}
                >
                  +{geneCount - 50} more
                </span>
              )}
            </div>
          ) : (
            <div style={{ fontSize: 12, color: "var(--text-secondary)" }}>
              No genes found
            </div>
          )}
        </div>
      )}

      {/* Top Neighbors */}
      {neighbors.length > 0 && (
        <div className="detail-section">
          <h4>Top Neighbors ({state.metric})</h4>
          {neighbors.map((n, i) => (
            <div key={i} className="detail-row">
              <span
                className="detail-label"
                style={{
                  maxWidth: 180,
                  overflow: "hidden",
                  textOverflow: "ellipsis",
                  whiteSpace: "nowrap",
                }}
              >
                {n.name}
              </span>
              <span className="detail-value">{n.weight.toFixed(1)}</span>
            </div>
          ))}
        </div>
      )}

      {/* Metric Profile */}
      <CollapsibleSection
        title="Metric Profile"
        loading={profileLoading}
        error={profileError}
        onExpand={fetchProfile}
        onRetry={fetchProfile}
        badge={metricProfile ? metricProfile.metrics.length : undefined}
      >
        {metricProfile && metricProfile.metrics.length > 0 && (
          <div style={{ display: "flex", flexDirection: "column", gap: 6 }}>
            {(() => {
              const maxAvg = Math.max(...metricProfile.metrics.map((m) => m.avg));
              return metricProfile.metrics.map((m) => (
                <div key={m.metric} style={{ fontSize: 11 }}>
                  <div
                    style={{
                      display: "flex",
                      justifyContent: "space-between",
                      marginBottom: 2,
                    }}
                  >
                    <span style={{ fontWeight: 500 }}>{m.metric}</span>
                    <span style={{ fontFamily: "var(--font-mono)", fontSize: 10 }}>
                      avg {m.avg.toFixed(2)} / max {m.max.toFixed(2)}
                    </span>
                  </div>
                  <div
                    style={{
                      height: 6,
                      background: "var(--bg-primary)",
                      borderRadius: 3,
                      overflow: "hidden",
                    }}
                  >
                    <div
                      style={{
                        height: "100%",
                        width: `${maxAvg > 0 ? (m.avg / maxAvg) * 100 : 0}%`,
                        background: "var(--accent)",
                        borderRadius: 3,
                      }}
                    />
                  </div>
                </div>
              ));
            })()}
          </div>
        )}
      </CollapsibleSection>

      {/* Hub Genes */}
      {state.info?.has_genes && (
        <CollapsibleSection
          title="Hub Genes"
          loading={hubLoading}
          error={hubError}
          onExpand={() => fetchHub(hubMethod)}
          onRetry={() => fetchHub(hubMethod)}
          badge={hubGenes ? hubGenes.length : undefined}
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
              ] as [HubMethod, string][]
            ).map(([method, label]) => (
              <button
                key={method}
                onClick={() => setHubMethod(method)}
                style={{
                  flex: 1,
                  padding: "4px 8px",
                  border: "none",
                  borderRadius: 3,
                  fontSize: 10,
                  fontWeight: 500,
                  cursor: "pointer",
                  background: hubMethod === method ? "var(--accent)" : "transparent",
                  color: hubMethod === method ? "white" : "var(--text-secondary)",
                }}
              >
                {label}
              </button>
            ))}
          </div>

          {/* Gene table */}
          {hubGenes && hubGenes.length > 0 && (
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
                  {hubGenes.map((g, i) => (
                    <tr
                      key={g.gene}
                      style={{
                        borderBottom:
                          i < hubGenes.length - 1
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

      {/* Neighborhood Enrichment */}
      {state.info?.has_genes && (
        <CollapsibleSection
          title="Neighborhood Enrichment"
          loading={enrichLoading}
          error={enrichError}
          onExpand={fetchEnrichment}
          onRetry={fetchEnrichment}
          badge={enrichGenes ? enrichGenes.length : undefined}
        >
          {enrichGenes && enrichGenes.length > 0 && (
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
                    <th style={{ padding: "5px 8px", textAlign: "right" }}>In/Total</th>
                  </tr>
                </thead>
                <tbody>
                  {enrichGenes.map((g, i) => (
                    <tr
                      key={g.gene}
                      style={{
                        borderBottom:
                          i < enrichGenes.length - 1
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
                      <td
                        style={{
                          padding: "5px 8px",
                          textAlign: "right",
                          fontSize: 10,
                          color: "var(--text-secondary)",
                        }}
                      >
                        {g.in_cluster_count}/{g.cluster_size}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}
        </CollapsibleSection>
      )}
    </>
  );
}
