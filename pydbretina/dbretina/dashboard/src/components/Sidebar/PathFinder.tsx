import { useState, useCallback, useEffect } from "react";
import { useDashboard } from "../../state/context";
import { fetchShortestPath, searchGroups } from "../../api";
import { useDebounce } from "../../hooks/useDebounce";

export default function PathFinder() {
  const { state, dispatch } = useDashboard();
  const [expanded, setExpanded] = useState(false);

  // Source search
  const [sourceQuery, setSourceQuery] = useState("");
  const [sourceResults, setSourceResults] = useState<{ id: number; name: string }[]>([]);
  const [selectedSource, setSelectedSource] = useState<string | null>(null);

  // Target search
  const [targetQuery, setTargetQuery] = useState("");
  const [targetResults, setTargetResults] = useState<{ id: number; name: string }[]>([]);
  const [selectedTarget, setSelectedTarget] = useState<string | null>(null);

  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [showGenes, setShowGenes] = useState(false);

  const debouncedSource = useDebounce(sourceQuery, 250);
  const debouncedTarget = useDebounce(targetQuery, 250);

  // Search for source groups
  const searchSource = useCallback(async () => {
    if (debouncedSource.length < 2) {
      setSourceResults([]);
      return;
    }
    try {
      const result = await searchGroups(debouncedSource, 10);
      setSourceResults(result.matches);
    } catch {
      setSourceResults([]);
    }
  }, [debouncedSource]);

  // Search for target groups
  const searchTarget = useCallback(async () => {
    if (debouncedTarget.length < 2) {
      setTargetResults([]);
      return;
    }
    try {
      const result = await searchGroups(debouncedTarget, 10);
      setTargetResults(result.matches);
    } catch {
      setTargetResults([]);
    }
  }, [debouncedTarget]);

  // Trigger searches on debounced query changes
  useEffect(() => {
    searchSource();
  }, [searchSource]);

  useEffect(() => {
    searchTarget();
  }, [searchTarget]);

  const handleFindPath = useCallback(async () => {
    if (!selectedSource || !selectedTarget) return;

    setLoading(true);
    setError(null);

    try {
      const result = await fetchShortestPath(
        selectedSource,
        selectedTarget,
        state.metric,
        state.cutoff
      );
      dispatch({ type: "SET_PATH_RESULT", result });

      // Highlight path nodes in graph
      if (result.connected && result.path_nodes.length > 0 && state.graphData) {
        const pathNames = new Set(result.path_nodes.map((n) => n.toLowerCase()));
        const nodeIds = state.graphData.nodes
          .filter((n) => pathNames.has(n.label.toLowerCase()))
          .map((n) => n.id);
        dispatch({ type: "SET_HIGHLIGHT_NODES", nodes: new Set(nodeIds) });
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : "Path finding failed");
      dispatch({ type: "SET_PATH_RESULT", result: null });
    } finally {
      setLoading(false);
    }
  }, [selectedSource, selectedTarget, state.metric, state.cutoff, state.graphData, dispatch]);

  const handleClear = useCallback(() => {
    dispatch({ type: "SET_PATH_RESULT", result: null });
    dispatch({ type: "SET_HIGHLIGHT_NODES", nodes: new Set() });
    setSelectedSource(null);
    setSelectedTarget(null);
    setSourceQuery("");
    setTargetQuery("");
  }, [dispatch]);

  const result = state.pathResult;

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
        <label style={{ cursor: "pointer", marginBottom: 0 }}>Path Finder</label>
        <span style={{ fontSize: 12, color: "var(--text-secondary)" }}>
          {expanded ? "[-]" : "[+]"}
        </span>
      </div>

      {expanded && (
        <div style={{ display: "flex", flexDirection: "column", gap: 12 }}>
          {/* Source input */}
          <div>
            <div
              style={{ fontSize: 11, color: "var(--text-secondary)", marginBottom: 4 }}
            >
              Source
            </div>
            {selectedSource ? (
              <div
                style={{
                  display: "flex",
                  alignItems: "center",
                  gap: 8,
                  padding: "6px 10px",
                  background: "var(--accent)",
                  color: "white",
                  borderRadius: 4,
                  fontSize: 12,
                }}
              >
                <span style={{ flex: 1 }}>{selectedSource}</span>
                <button
                  onClick={() => {
                    setSelectedSource(null);
                    setSourceQuery("");
                  }}
                  style={{
                    background: "none",
                    border: "none",
                    color: "white",
                    cursor: "pointer",
                    fontSize: 14,
                  }}
                >
                  ×
                </button>
              </div>
            ) : (
              <div style={{ position: "relative" }}>
                <input
                  type="text"
                  value={sourceQuery}
                  onChange={(e) => {
                    setSourceQuery(e.target.value);
                    if (e.target.value.length >= 2) searchSource();
                  }}
                  placeholder="Search source group..."
                  style={{
                    width: "100%",
                    padding: "6px 10px",
                    border: "1px solid var(--border-color)",
                    borderRadius: 4,
                    fontSize: 12,
                  }}
                />
                {sourceResults.length > 0 && (
                  <div
                    style={{
                      position: "absolute",
                      top: "100%",
                      left: 0,
                      right: 0,
                      background: "white",
                      border: "1px solid var(--border-color)",
                      borderRadius: 4,
                      maxHeight: 150,
                      overflowY: "auto",
                      zIndex: 10,
                    }}
                  >
                    {sourceResults.map((r) => (
                      <div
                        key={r.id}
                        onClick={() => {
                          setSelectedSource(r.name);
                          setSourceQuery("");
                          setSourceResults([]);
                        }}
                        style={{
                          padding: "6px 10px",
                          fontSize: 12,
                          cursor: "pointer",
                        }}
                        onMouseOver={(e) =>
                          (e.currentTarget.style.background = "var(--bg-primary)")
                        }
                        onMouseOut={(e) => (e.currentTarget.style.background = "white")}
                      >
                        {r.name}
                      </div>
                    ))}
                  </div>
                )}
              </div>
            )}
          </div>

          {/* Target input */}
          <div>
            <div
              style={{ fontSize: 11, color: "var(--text-secondary)", marginBottom: 4 }}
            >
              Target
            </div>
            {selectedTarget ? (
              <div
                style={{
                  display: "flex",
                  alignItems: "center",
                  gap: 8,
                  padding: "6px 10px",
                  background: "var(--accent)",
                  color: "white",
                  borderRadius: 4,
                  fontSize: 12,
                }}
              >
                <span style={{ flex: 1 }}>{selectedTarget}</span>
                <button
                  onClick={() => {
                    setSelectedTarget(null);
                    setTargetQuery("");
                  }}
                  style={{
                    background: "none",
                    border: "none",
                    color: "white",
                    cursor: "pointer",
                    fontSize: 14,
                  }}
                >
                  ×
                </button>
              </div>
            ) : (
              <div style={{ position: "relative" }}>
                <input
                  type="text"
                  value={targetQuery}
                  onChange={(e) => {
                    setTargetQuery(e.target.value);
                    if (e.target.value.length >= 2) searchTarget();
                  }}
                  placeholder="Search target group..."
                  style={{
                    width: "100%",
                    padding: "6px 10px",
                    border: "1px solid var(--border-color)",
                    borderRadius: 4,
                    fontSize: 12,
                  }}
                />
                {targetResults.length > 0 && (
                  <div
                    style={{
                      position: "absolute",
                      top: "100%",
                      left: 0,
                      right: 0,
                      background: "white",
                      border: "1px solid var(--border-color)",
                      borderRadius: 4,
                      maxHeight: 150,
                      overflowY: "auto",
                      zIndex: 10,
                    }}
                  >
                    {targetResults.map((r) => (
                      <div
                        key={r.id}
                        onClick={() => {
                          setSelectedTarget(r.name);
                          setTargetQuery("");
                          setTargetResults([]);
                        }}
                        style={{
                          padding: "6px 10px",
                          fontSize: 12,
                          cursor: "pointer",
                        }}
                        onMouseOver={(e) =>
                          (e.currentTarget.style.background = "var(--bg-primary)")
                        }
                        onMouseOut={(e) => (e.currentTarget.style.background = "white")}
                      >
                        {r.name}
                      </div>
                    ))}
                  </div>
                )}
              </div>
            )}
          </div>

          {/* Find path button */}
          <div style={{ display: "flex", gap: 8 }}>
            <button
              onClick={handleFindPath}
              disabled={!selectedSource || !selectedTarget || loading}
              style={{
                flex: 1,
                padding: "10px 16px",
                background:
                  selectedSource && selectedTarget && !loading ? "var(--accent)" : "#ccc",
                color: "white",
                border: "none",
                borderRadius: 4,
                fontSize: 13,
                fontWeight: 500,
                cursor:
                  selectedSource && selectedTarget && !loading ? "pointer" : "not-allowed",
              }}
            >
              {loading ? "Finding..." : "Find Path"}
            </button>
            {result && (
              <button
                onClick={handleClear}
                style={{
                  padding: "10px 12px",
                  background: "white",
                  color: "var(--text-primary)",
                  border: "1px solid var(--border-color)",
                  borderRadius: 4,
                  fontSize: 13,
                  cursor: "pointer",
                }}
              >
                Clear
              </button>
            )}
          </div>

          {error && <div style={{ fontSize: 12, color: "#e74c3c" }}>{error}</div>}

          {/* Path result */}
          {result && (
            <div
              style={{
                padding: 12,
                background: "var(--bg-primary)",
                borderRadius: 6,
                border: "1px solid var(--border-color)",
              }}
            >
              {result.connected ? (
                <>
                  <div
                    style={{
                      fontSize: 14,
                      fontWeight: 600,
                      marginBottom: 8,
                      color: "var(--accent)",
                    }}
                  >
                    Path Found! ({result.path_length} hops)
                  </div>

                  {/* Path visualization */}
                  <div style={{ display: "flex", flexDirection: "column", gap: 4 }}>
                    {result.path_nodes.map((node, i) => (
                      <div key={i}>
                        <div
                          style={{
                            padding: "6px 10px",
                            background: "white",
                            borderRadius: 4,
                            fontSize: 12,
                            fontWeight: i === 0 || i === result.path_nodes.length - 1 ? 600 : 400,
                            border:
                              i === 0 || i === result.path_nodes.length - 1
                                ? "2px solid var(--accent)"
                                : "1px solid var(--border-color)",
                          }}
                        >
                          {node}
                        </div>
                        {i < result.path_nodes.length - 1 && (
                          <div
                            style={{
                              textAlign: "center",
                              fontSize: 16,
                              color: "var(--text-secondary)",
                              padding: "2px 0",
                            }}
                          >
                            ↓
                          </div>
                        )}
                      </div>
                    ))}
                  </div>

                  {/* Shared genes along path */}
                  {result.shared_genes_along_path &&
                    result.shared_genes_along_path.length > 0 && (
                      <div style={{ marginTop: 12 }}>
                        <button
                          onClick={() => setShowGenes(!showGenes)}
                          style={{
                            width: "100%",
                            padding: "6px 10px",
                            background: "white",
                            border: "1px solid var(--border-color)",
                            borderRadius: 4,
                            fontSize: 11,
                            cursor: "pointer",
                            textAlign: "left",
                          }}
                        >
                          {showGenes ? "Hide" : "Show"} Shared Genes Along Path
                        </button>

                        {showGenes && (
                          <div style={{ marginTop: 8 }}>
                            {result.shared_genes_along_path.map((edge, i) => (
                              <div
                                key={i}
                                style={{
                                  padding: 8,
                                  background: "white",
                                  borderRadius: 4,
                                  marginTop: i > 0 ? 8 : 0,
                                  fontSize: 11,
                                }}
                              >
                                <div style={{ fontWeight: 600, marginBottom: 4 }}>
                                  {edge.from} → {edge.to} ({edge.shared_count} genes)
                                </div>
                                <div
                                  style={{
                                    display: "flex",
                                    flexWrap: "wrap",
                                    gap: 4,
                                  }}
                                >
                                  {edge.genes.slice(0, 10).map((gene) => (
                                    <span
                                      key={gene}
                                      style={{
                                        padding: "2px 6px",
                                        background: "var(--bg-primary)",
                                        borderRadius: 3,
                                        fontSize: 10,
                                      }}
                                    >
                                      {gene}
                                    </span>
                                  ))}
                                  {edge.genes.length > 10 && (
                                    <span
                                      style={{
                                        padding: "2px 6px",
                                        color: "var(--text-secondary)",
                                        fontSize: 10,
                                      }}
                                    >
                                      +{edge.genes.length - 10} more
                                    </span>
                                  )}
                                </div>
                              </div>
                            ))}
                          </div>
                        )}
                      </div>
                    )}
                </>
              ) : (
                <div style={{ fontSize: 13, color: "#e74c3c" }}>
                  No path found between these groups at the current cutoff.
                  <div style={{ fontSize: 11, marginTop: 4, color: "var(--text-secondary)" }}>
                    Try lowering the cutoff value to include weaker connections.
                  </div>
                </div>
              )}
            </div>
          )}

          {/* Help text */}
          {!result && !loading && (
            <div style={{ fontSize: 11, color: "var(--text-secondary)" }}>
              Find the shortest path between two diseases/groups in the network. See which
              intermediate diseases connect them and what genes they share.
            </div>
          )}
        </div>
      )}
    </div>
  );
}
