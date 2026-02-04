import { useState, useEffect, useCallback } from "react";
import { useDashboard } from "../../state/context";
import { fetchSharedFeatures, explainPair, ExplainPairResult } from "../../api";

export default function EdgeDetail() {
  const { state, dispatch } = useDashboard();
  const edge = state.selectedEdge;
  const [genes, setGenes] = useState<string[]>([]);
  const [geneCount, setGeneCount] = useState<number>(0);
  const [loading, setLoading] = useState(false);

  // Explain pair state
  const [showExplain, setShowExplain] = useState(false);
  const [explainLoading, setExplainLoading] = useState(false);
  const [explainResult, setExplainResult] = useState<ExplainPairResult | null>(null);
  const [explainError, setExplainError] = useState<string | null>(null);

  if (!edge) return null;

  const srcNode = state.graphData?.nodes.find((n) => n.id === edge.source);
  const dstNode = state.graphData?.nodes.find((n) => n.id === edge.target);
  const edgeData = state.graphData?.edges.find(
    (e) =>
      (e.source === edge.source && e.target === edge.target) ||
      (e.source === edge.target && e.target === edge.source)
  );

  const srcName = srcNode?.label ?? edge.source;
  const dstName = dstNode?.label ?? edge.target;

  // Fetch shared genes
  useEffect(() => {
    if (!state.info?.has_genes || !srcNode || !dstNode) return;
    setLoading(true);
    fetchSharedFeatures(srcName, dstName)
      .then((res) => {
        setGenes(res.features.slice(0, 50));
        setGeneCount(res.count);
      })
      .catch(() => setGenes([]))
      .finally(() => setLoading(false));
  }, [srcName, dstName, state.info?.has_genes]);

  // Reset explain state when edge changes
  useEffect(() => {
    setShowExplain(false);
    setExplainResult(null);
    setExplainError(null);
  }, [srcName, dstName]);

  const handleExplainConnection = useCallback(async () => {
    if (!srcName || !dstName) return;

    setShowExplain(true);
    setExplainLoading(true);
    setExplainError(null);

    try {
      const result = await explainPair(srcName, dstName, "hypergraph");
      setExplainResult(result);
      dispatch({ type: "SET_EXPLAIN_PAIR", result });
    } catch (err) {
      setExplainError(err instanceof Error ? err.message : "Analysis failed");
      setExplainResult(null);
    } finally {
      setExplainLoading(false);
    }
  }, [srcName, dstName, dispatch]);

  return (
    <>
      <div className="detail-section">
        <h4>Connected Groups</h4>
        <div style={{ fontSize: 12, marginBottom: 4 }}>
          <strong>{srcName}</strong>
        </div>
        <div style={{ fontSize: 11, color: "var(--text-secondary)", marginBottom: 4 }}>
          &harr;
        </div>
        <div style={{ fontSize: 12 }}>
          <strong>{dstName}</strong>
        </div>
      </div>
      {edgeData && (
        <div className="detail-section">
          <h4>Metrics</h4>
          <div className="detail-row">
            <span className="detail-label">Weight ({state.metric})</span>
            <span className="detail-value">{edgeData.weight.toFixed(2)}</span>
          </div>
          <div className="detail-row">
            <span className="detail-label">Shared features</span>
            <span className="detail-value">{edgeData.shared_features}</span>
          </div>
        </div>
      )}
      {state.info?.has_genes && (
        <div className="detail-section">
          <h4>Shared Genes {geneCount > 0 && `(${geneCount})`}</h4>
          {loading ? (
            <div style={{ fontSize: 12, color: "var(--text-secondary)" }}>Loading...</div>
          ) : genes.length > 0 ? (
            <div className="gene-list">
              {genes.map((g) => (
                <span key={g} className="gene-tag">
                  {g}
                </span>
              ))}
              {geneCount > 50 && (
                <span className="gene-tag" style={{ color: "var(--text-secondary)" }}>
                  +{geneCount - 50} more
                </span>
              )}
            </div>
          ) : (
            <div style={{ fontSize: 12, color: "var(--text-secondary)" }}>No shared genes found</div>
          )}
        </div>
      )}

      {/* Why Connected? - Gene importance analysis */}
      {state.info?.has_genes && geneCount > 0 && (
        <div className="detail-section">
          {!showExplain ? (
            <button
              onClick={handleExplainConnection}
              style={{
                width: "100%",
                padding: "10px 16px",
                background: "var(--accent)",
                color: "white",
                border: "none",
                borderRadius: 4,
                fontSize: 13,
                fontWeight: 500,
                cursor: "pointer",
              }}
            >
              Why Connected?
            </button>
          ) : (
            <div>
              <div
                style={{
                  display: "flex",
                  justifyContent: "space-between",
                  alignItems: "center",
                  marginBottom: 8,
                }}
              >
                <h4 style={{ margin: 0 }}>Important Genes</h4>
                <button
                  onClick={() => {
                    setShowExplain(false);
                    setExplainResult(null);
                  }}
                  style={{
                    padding: "2px 8px",
                    fontSize: 10,
                    border: "1px solid var(--border-color)",
                    borderRadius: 3,
                    background: "white",
                    cursor: "pointer",
                  }}
                >
                  Hide
                </button>
              </div>

              {explainLoading ? (
                <div style={{ fontSize: 12, color: "var(--text-secondary)" }}>
                  Analyzing connection...
                </div>
              ) : explainError ? (
                <div style={{ fontSize: 12, color: "#e74c3c" }}>{explainError}</div>
              ) : explainResult ? (
                <div>
                  <div style={{ fontSize: 10, color: "var(--text-secondary)", marginBottom: 8 }}>
                    Genes ranked by specificity to this connection (TF-IDF method)
                  </div>
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
                          <th style={{ padding: "6px 8px", textAlign: "left" }}>Gene</th>
                          <th style={{ padding: "6px 8px", textAlign: "right" }}>Score</th>
                        </tr>
                      </thead>
                      <tbody>
                        {explainResult.genes.slice(0, 30).map((gene, i) => (
                          <tr
                            key={gene.gene}
                            style={{
                              borderBottom:
                                i < Math.min(explainResult.genes.length, 30) - 1
                                  ? "1px solid var(--border-color)"
                                  : "none",
                            }}
                          >
                            <td
                              style={{
                                padding: "6px 8px",
                                fontWeight: i < 5 ? 600 : 400,
                                color: i < 5 ? "var(--accent)" : "inherit",
                              }}
                            >
                              {gene.gene}
                              {gene.specificity !== undefined && gene.specificity > 0.8 && (
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
                            <td style={{ padding: "6px 8px", textAlign: "right" }}>
                              {gene.score.toFixed(4)}
                            </td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                  {explainResult.genes.length > 30 && (
                    <div
                      style={{
                        fontSize: 10,
                        color: "var(--text-secondary)",
                        marginTop: 4,
                        textAlign: "center",
                      }}
                    >
                      Showing top 30 of {explainResult.genes.length} genes
                    </div>
                  )}
                </div>
              ) : null}
            </div>
          )}
        </div>
      )}
    </>
  );
}
