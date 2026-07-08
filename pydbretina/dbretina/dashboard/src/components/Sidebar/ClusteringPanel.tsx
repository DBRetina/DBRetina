import { useState, useEffect, useCallback } from "react";
import { useDashboard } from "../../state/context";
import {
  fetchClusteringAlgorithms,
  runClustering,
  AlgorithmInfo,
  AlgorithmParameter,
  ClusteringResult,
  getAPIError,
} from "../../api";
import { LoadingState, EmptyState } from "../Common";

interface ParameterInputProps {
  param: AlgorithmParameter;
  value: unknown;
  onChange: (value: unknown) => void;
}

function ParameterInput({ param, value, onChange }: ParameterInputProps) {
  const currentValue = value ?? param.default;

  if (param.type === "bool") {
    return (
      <label style={{ display: "flex", alignItems: "center", gap: 8, cursor: "pointer" }}>
        <input
          type="checkbox"
          checked={Boolean(currentValue)}
          onChange={(e) => onChange(e.target.checked)}
          style={{ width: 16, height: 16 }}
        />
        <span style={{ fontSize: 12 }}>{param.description}</span>
      </label>
    );
  }

  if (param.choices && param.choices.length > 0) {
    return (
      <select
        value={String(currentValue)}
        onChange={(e) => onChange(e.target.value)}
        style={{
          width: "100%",
          padding: "6px 10px",
          border: "1px solid var(--border-color)",
          borderRadius: 4,
          fontSize: 13,
          background: "var(--bg-primary)",
        }}
      >
        {param.choices.map((choice) => (
          <option key={choice} value={choice}>
            {choice}
          </option>
        ))}
      </select>
    );
  }

  if (param.type === "int" || param.type === "float") {
    const step = param.type === "float" ? 0.1 : 1;
    return (
      <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
        <input
          type="range"
          min={param.min_value ?? 0}
          max={param.max_value ?? (param.type === "int" ? 100 : 10)}
          step={step}
          value={Number(currentValue)}
          onChange={(e) => onChange(param.type === "int" ? parseInt(e.target.value) : parseFloat(e.target.value))}
          style={{ flex: 1, accentColor: "var(--accent)" }}
        />
        <span style={{ minWidth: 40, textAlign: "right", fontSize: 12, fontFamily: "var(--font-mono)" }}>
          {param.type === "float" ? Number(currentValue).toFixed(2) : String(currentValue)}
        </span>
      </div>
    );
  }

  return (
    <input
      type="text"
      value={String(currentValue)}
      onChange={(e) => onChange(e.target.value)}
      style={{
        width: "100%",
        padding: "6px 10px",
        border: "1px solid var(--border-color)",
        borderRadius: 4,
        fontSize: 13,
        background: "var(--bg-primary)",
      }}
    />
  );
}

export default function ClusteringPanel() {
  const { state, dispatch } = useDashboard();

  // Algorithm state
  const [algorithms, setAlgorithms] = useState<AlgorithmInfo[]>([]);
  const [selectedAlgorithm, setSelectedAlgorithm] = useState<string>("leiden");
  const [parameters, setParameters] = useState<Record<string, unknown>>({});

  // UI state
  const [loading, setLoading] = useState(false);
  const [algorithmsLoading, setAlgorithmsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<ClusteringResult | null>(null);
  const [expanded, setExpanded] = useState(false);

  // Load algorithms on mount
  useEffect(() => {
    fetchClusteringAlgorithms()
      .then((data) => {
        setAlgorithms(data.algorithms);
        // Set default parameters for initial algorithm
        const defaultAlgo = data.algorithms.find((a) => a.name === "leiden");
        if (defaultAlgo) {
          const defaultParams: Record<string, unknown> = {};
          defaultAlgo.parameters.forEach((p) => {
            defaultParams[p.name] = p.default;
          });
          setParameters(defaultParams);
        }
      })
      .catch((e) => setError(getAPIError(e)?.detail || e.message))
      .finally(() => setAlgorithmsLoading(false));
  }, []);

  // Update parameters when algorithm changes
  const handleAlgorithmChange = useCallback(
    (name: string) => {
      setSelectedAlgorithm(name);
      const algo = algorithms.find((a) => a.name === name);
      if (algo) {
        const defaultParams: Record<string, unknown> = {};
        algo.parameters.forEach((p) => {
          defaultParams[p.name] = p.default;
        });
        setParameters(defaultParams);
      }
      setResult(null);
      setError(null);
    },
    [algorithms]
  );

  // Run clustering
  const handleRunClustering = useCallback(async () => {
    setLoading(true);
    setError(null);

    try {
      const clusterResult = await runClustering(
        selectedAlgorithm,
        parameters,
        state.metric,
        state.cutoff
      );
      setResult(clusterResult);

      // Update graph highlighting based on clusters (optional)
      // dispatch({ type: "SET_CLUSTER_RESULT", result: clusterResult });
    } catch (e) {
      const apiError = getAPIError(e);
      setError(apiError?.detail || (e instanceof Error ? e.message : "Clustering failed"));
    } finally {
      setLoading(false);
    }
  }, [selectedAlgorithm, parameters, state.metric, state.cutoff]);

  const selectedAlgoInfo = algorithms.find((a) => a.name === selectedAlgorithm);

  if (algorithmsLoading) {
    return (
      <div className="sidebar-section">
        <label>Clustering</label>
        <LoadingState message="Loading algorithms..." size="small" />
      </div>
    );
  }

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
        <label style={{ cursor: "pointer", marginBottom: 0 }}>Clustering</label>
        <span style={{ fontSize: 12, color: "var(--text-secondary)" }}>
          {expanded ? "[-]" : "[+]"}
        </span>
      </div>

      {expanded && (
        <div style={{ display: "flex", flexDirection: "column", gap: 12 }}>
          {/* Algorithm selector */}
          <div>
            <div style={{ fontSize: 11, color: "var(--text-secondary)", marginBottom: 4 }}>
              Algorithm
            </div>
            <select
              value={selectedAlgorithm}
              onChange={(e) => handleAlgorithmChange(e.target.value)}
              style={{
                width: "100%",
                padding: "6px 10px",
                border: "1px solid var(--border-color)",
                borderRadius: 4,
                fontSize: 13,
                background: "var(--bg-primary)",
              }}
            >
              {algorithms.map((algo) => (
                <option key={algo.name} value={algo.name}>
                  {algo.display_name}
                </option>
              ))}
            </select>
          </div>

          {/* Algorithm description */}
          {selectedAlgoInfo && (
            <div style={{ fontSize: 11, color: "var(--text-secondary)", lineHeight: 1.4 }}>
              {selectedAlgoInfo.description}
            </div>
          )}

          {/* Parameters */}
          {selectedAlgoInfo && selectedAlgoInfo.parameters.length > 0 && (
            <div style={{ display: "flex", flexDirection: "column", gap: 10 }}>
              <div style={{ fontSize: 11, color: "var(--text-secondary)", fontWeight: 600 }}>
                Parameters
              </div>
              {selectedAlgoInfo.parameters.map((param) => (
                <div key={param.name}>
                  <div
                    style={{ fontSize: 11, color: "var(--text-secondary)", marginBottom: 4 }}
                    title={param.description}
                  >
                    {param.name.replace(/_/g, " ")}
                  </div>
                  <ParameterInput
                    param={param}
                    value={parameters[param.name]}
                    onChange={(value) =>
                      setParameters((prev) => ({ ...prev, [param.name]: value }))
                    }
                  />
                </div>
              ))}
            </div>
          )}

          {/* Run button */}
          <button
            onClick={handleRunClustering}
            disabled={loading || !state.graphData}
            style={{
              padding: "8px 16px",
              background: loading ? "var(--text-secondary)" : "var(--accent)",
              color: "white",
              border: "none",
              borderRadius: 4,
              fontSize: 12,
              fontWeight: 500,
              cursor: loading ? "wait" : "pointer",
            }}
          >
            {loading ? "Running..." : "Run Clustering"}
          </button>

          {/* Error display */}
          {error && (
            <div
              style={{
                padding: 8,
                background: "#ffeaed",
                color: "var(--danger)",
                borderRadius: 4,
                fontSize: 11,
              }}
            >
              {error}
            </div>
          )}

          {/* Results */}
          {result && (
            <div
              style={{
                padding: 12,
                background: "var(--bg-primary)",
                borderRadius: 4,
                fontSize: 12,
              }}
            >
              <div style={{ fontWeight: 600, marginBottom: 8 }}>Results</div>

              <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 8 }}>
                <div>
                  <div style={{ color: "var(--text-secondary)", fontSize: 10 }}>Clusters</div>
                  <div style={{ fontWeight: 600, fontFamily: "var(--font-mono)" }}>
                    {result.num_clusters}
                  </div>
                </div>
                <div>
                  <div style={{ color: "var(--text-secondary)", fontSize: 10 }}>Modularity</div>
                  <div style={{ fontWeight: 600, fontFamily: "var(--font-mono)" }}>
                    {result.modularity !== null ? result.modularity.toFixed(4) : "N/A"}
                  </div>
                </div>
              </div>

              {/* Cluster sizes distribution */}
              <div style={{ marginTop: 12 }}>
                <div style={{ color: "var(--text-secondary)", fontSize: 10, marginBottom: 4 }}>
                  Cluster Sizes (top 5)
                </div>
                <div style={{ display: "flex", flexDirection: "column", gap: 4 }}>
                  {Object.entries(result.cluster_sizes)
                    .sort(([, a], [, b]) => b - a)
                    .slice(0, 5)
                    .map(([clusterId, size]) => (
                      <div
                        key={clusterId}
                        style={{
                          display: "flex",
                          alignItems: "center",
                          gap: 8,
                        }}
                      >
                        <div
                          style={{
                            width: 10,
                            height: 10,
                            borderRadius: "50%",
                            background: `hsl(${(parseInt(clusterId) * 137.5) % 360}, 70%, 50%)`,
                          }}
                        />
                        <div style={{ flex: 1, fontSize: 11 }}>Cluster {clusterId}</div>
                        <div style={{ fontFamily: "var(--font-mono)", fontSize: 11 }}>{size}</div>
                      </div>
                    ))}
                </div>
              </div>

              {/* Export button */}
              <button
                onClick={() => {
                  const blob = new Blob([JSON.stringify(result, null, 2)], {
                    type: "application/json",
                  });
                  const url = URL.createObjectURL(blob);
                  const a = document.createElement("a");
                  a.href = url;
                  a.download = `clustering_${result.algorithm}_${Date.now()}.json`;
                  a.click();
                  URL.revokeObjectURL(url);
                }}
                style={{
                  marginTop: 12,
                  padding: "6px 12px",
                  background: "transparent",
                  color: "var(--accent)",
                  border: "1px solid var(--accent)",
                  borderRadius: 4,
                  fontSize: 11,
                  cursor: "pointer",
                  width: "100%",
                }}
              >
                Export Results (JSON)
              </button>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
