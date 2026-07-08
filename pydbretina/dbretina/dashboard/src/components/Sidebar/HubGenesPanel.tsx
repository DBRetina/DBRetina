import { useState, useCallback } from "react";
import { useDashboard } from "../../state/context";
import { fetchHubGenes } from "../../api";

type ScoringMethod = "hypergraph" | "edge_weighted" | "projection";

const METHODS: { value: ScoringMethod; label: string; description: string }[] = [
  {
    value: "hypergraph",
    label: "TF-IDF",
    description: "Genes specific to this disease neighborhood (recommended)",
  },
  {
    value: "edge_weighted",
    label: "Edge Weighted",
    description: "Genes weighted by similarity strength",
  },
  {
    value: "projection",
    label: "PageRank",
    description: "Central genes in gene co-occurrence network",
  },
];

export default function HubGenesPanel() {
  const { state, dispatch } = useDashboard();
  const [expanded, setExpanded] = useState(false);
  const [method, setMethod] = useState<ScoringMethod>("hypergraph");
  const [hops, setHops] = useState(2);
  const [topN, setTopN] = useState(30);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const selectedGroup = state.selectedNode?.label || state.focusGroup;

  const handleAnalyze = useCallback(async () => {
    if (!selectedGroup) return;

    setLoading(true);
    setError(null);

    try {
      const result = await fetchHubGenes(selectedGroup, method, hops, topN);
      dispatch({ type: "SET_HUB_GENES", result });
    } catch (err) {
      setError(err instanceof Error ? err.message : "Analysis failed");
      dispatch({ type: "SET_HUB_GENES", result: null });
    } finally {
      setLoading(false);
    }
  }, [selectedGroup, method, hops, topN, dispatch]);

  const handleClear = useCallback(() => {
    dispatch({ type: "SET_HUB_GENES", result: null });
  }, [dispatch]);

  const result = state.hubGenesResult;

  return (
    <div className="sidebar-section">
      <div
        onClick={() => setExpanded(!expanded)}
        style={{
          display: "flex",
          alignItems: "center",
          justifyContent: "space-between",
          cursor: "pointer",
          marginBottom: expanded ? 12 : 0,
        }}
      >
        <label style={{ cursor: "pointer", marginBottom: 0 }}>Hub Genes</label>
        <span style={{ fontSize: 12, color: "var(--text-secondary)" }}>
          {expanded ? "[-]" : "[+]"}
        </span>
      </div>

      {expanded && (
        <div style={{ display: "flex", flexDirection: "column", gap: 12 }}>
          {/* Target group info */}
          {selectedGroup ? (
            <div
              style={{
                padding: 8,
                background: "var(--bg-primary)",
                borderRadius: 4,
                fontSize: 12,
              }}
            >
              <span style={{ color: "var(--text-secondary)" }}>Analyzing: </span>
              <strong>{selectedGroup}</strong>
            </div>
          ) : (
            <div
              style={{
                padding: 8,
                background: "#fff3cd",
                borderRadius: 4,
                fontSize: 12,
                color: "#856404",
              }}
            >
              Select a node or focus on a group to analyze hub genes
            </div>
          )}

          {/* Method selector */}
          <div>
            <div
              style={{
                fontSize: 11,
                color: "var(--text-secondary)",
                marginBottom: 4,
                fontWeight: 500,
              }}
            >
              Scoring Method
            </div>
            <select
              value={method}
              onChange={(e) => setMethod(e.target.value as ScoringMethod)}
              style={{
                width: "100%",
                padding: "6px 8px",
                fontSize: 12,
                border: "1px solid var(--border-color)",
                borderRadius: 4,
              }}
            >
              {METHODS.map((m) => (
                <option key={m.value} value={m.value}>
                  {m.label}
                </option>
              ))}
            </select>
            <div style={{ fontSize: 10, color: "var(--text-secondary)", marginTop: 2 }}>
              {METHODS.find((m) => m.value === method)?.description}
            </div>
          </div>

          {/* Parameters */}
          <div style={{ display: "flex", gap: 12 }}>
            <div style={{ flex: 1 }}>
              <div
                style={{
                  fontSize: 11,
                  color: "var(--text-secondary)",
                  marginBottom: 4,
                }}
              >
                Hops: {hops}
              </div>
              <input
                type="range"
                min={1}
                max={5}
                value={hops}
                onChange={(e) => setHops(parseInt(e.target.value))}
                style={{ width: "100%" }}
              />
            </div>
            <div style={{ flex: 1 }}>
              <div
                style={{
                  fontSize: 11,
                  color: "var(--text-secondary)",
                  marginBottom: 4,
                }}
              >
                Top N: {topN}
              </div>
              <input
                type="range"
                min={10}
                max={100}
                step={10}
                value={topN}
                onChange={(e) => setTopN(parseInt(e.target.value))}
                style={{ width: "100%" }}
              />
            </div>
          </div>

          {/* Analyze button */}
          <button
            onClick={handleAnalyze}
            disabled={!selectedGroup || loading}
            style={{
              padding: "10px 16px",
              background: selectedGroup && !loading ? "var(--accent)" : "#ccc",
              color: "white",
              border: "none",
              borderRadius: 4,
              fontSize: 13,
              fontWeight: 500,
              cursor: selectedGroup && !loading ? "pointer" : "not-allowed",
            }}
          >
            {loading ? "Analyzing..." : "Find Hub Genes"}
          </button>

          {error && (
            <div style={{ fontSize: 12, color: "#e74c3c" }}>{error}</div>
          )}

          {/* Results */}
          {result && (
            <div>
              <div
                style={{
                  display: "flex",
                  justifyContent: "space-between",
                  alignItems: "center",
                  marginBottom: 8,
                }}
              >
                <div
                  style={{
                    fontSize: 11,
                    fontWeight: 600,
                    color: "var(--text-secondary)",
                    textTransform: "uppercase",
                  }}
                >
                  Hub Genes for {result.group}
                </div>
                <button
                  onClick={handleClear}
                  style={{
                    padding: "2px 8px",
                    fontSize: 10,
                    border: "1px solid var(--border-color)",
                    borderRadius: 3,
                    background: "white",
                    cursor: "pointer",
                  }}
                >
                  Clear
                </button>
              </div>

              <div
                style={{
                  maxHeight: 250,
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
                      <th style={{ padding: "6px 8px", textAlign: "right" }}>Diseases</th>
                    </tr>
                  </thead>
                  <tbody>
                    {result.genes.map((gene, i) => (
                      <tr
                        key={gene.gene}
                        style={{
                          borderBottom:
                            i < result.genes.length - 1 ? "1px solid var(--border-color)" : "none",
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
                        </td>
                        <td style={{ padding: "6px 8px", textAlign: "right" }}>
                          {gene.score.toFixed(4)}
                        </td>
                        <td style={{ padding: "6px 8px", textAlign: "right" }}>
                          {gene.num_diseases || "-"}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>

              <div style={{ fontSize: 10, color: "var(--text-secondary)", marginTop: 4 }}>
                Method: {result.method}, Hops: {result.hops}
              </div>
            </div>
          )}

          {/* Help text */}
          {!result && !loading && (
            <div style={{ fontSize: 11, color: "var(--text-secondary)" }}>
              Identify the most important genes driving connections for a disease. Hub genes are
              ranked by their importance in the disease neighborhood.
            </div>
          )}
        </div>
      )}
    </div>
  );
}
