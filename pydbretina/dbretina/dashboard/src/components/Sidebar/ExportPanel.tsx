import { useState, useCallback } from "react";
import { useDashboard } from "../../state/context";
import {
  getDataExportUrl,
  getGraphExportUrl,
  DataExportFormat,
  GraphExportFormat,
} from "../../api";

export default function ExportPanel() {
  const { state } = useDashboard();
  const [expanded, setExpanded] = useState(false);

  // Data export options
  const [dataFormat, setDataFormat] = useState<DataExportFormat>("csv");
  const [dataLimit, setDataLimit] = useState(100000);

  // Graph export options
  const [graphFormat, setGraphFormat] = useState<GraphExportFormat>("json");
  const [includeLayout, setIncludeLayout] = useState(false);
  const [includeCommunities, setIncludeCommunities] = useState(true);

  const handleDataExport = useCallback(() => {
    if (!state.info) return;
    const url = getDataExportUrl(dataFormat, state.metric, state.cutoff, dataLimit);
    window.open(url, "_blank");
  }, [dataFormat, state.metric, state.cutoff, dataLimit, state.info]);

  const handleGraphExport = useCallback(() => {
    if (!state.info) return;
    const url = getGraphExportUrl(
      graphFormat,
      state.metric,
      state.cutoff,
      includeLayout,
      includeCommunities
    );
    window.open(url, "_blank");
  }, [graphFormat, state.metric, state.cutoff, includeLayout, includeCommunities, state.info]);

  const dataFormats: { value: DataExportFormat; label: string; desc: string }[] = [
    { value: "csv", label: "CSV", desc: "Comma-separated values" },
    { value: "tsv", label: "TSV", desc: "Tab-separated values" },
    { value: "json", label: "JSON", desc: "JavaScript Object Notation" },
  ];

  const graphFormats: { value: GraphExportFormat; label: string; desc: string }[] = [
    { value: "json", label: "Cytoscape JSON", desc: "For Cytoscape.js visualization" },
    { value: "graphml", label: "GraphML", desc: "Standard graph format" },
    { value: "gexf", label: "GEXF", desc: "Graph Exchange XML Format" },
  ];

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
        <label style={{ cursor: "pointer", marginBottom: 0 }}>Export</label>
        <span style={{ fontSize: 12, color: "var(--text-secondary)" }}>
          {expanded ? "[-]" : "[+]"}
        </span>
      </div>

      {expanded && (
        <div style={{ display: "flex", flexDirection: "column", gap: 16 }}>
          {/* Data Export Section */}
          <div>
            <div
              style={{
                fontSize: 11,
                fontWeight: 600,
                color: "var(--text-secondary)",
                marginBottom: 8,
                textTransform: "uppercase",
                letterSpacing: 0.5,
              }}
            >
              Pairwise Data
            </div>

            {/* Format selector */}
            <div style={{ display: "flex", gap: 4, marginBottom: 8 }}>
              {dataFormats.map((fmt) => (
                <button
                  key={fmt.value}
                  onClick={() => setDataFormat(fmt.value)}
                  title={fmt.desc}
                  style={{
                    flex: 1,
                    padding: "6px 8px",
                    fontSize: 11,
                    border: `1px solid ${dataFormat === fmt.value ? "var(--accent)" : "var(--border-color)"}`,
                    borderRadius: 4,
                    background: dataFormat === fmt.value ? "var(--accent)" : "white",
                    color: dataFormat === fmt.value ? "white" : "var(--text-primary)",
                    cursor: "pointer",
                  }}
                >
                  {fmt.label}
                </button>
              ))}
            </div>

            {/* Limit input */}
            <div style={{ marginBottom: 8 }}>
              <div style={{ fontSize: 10, color: "var(--text-secondary)", marginBottom: 4 }}>
                Max rows
              </div>
              <input
                type="number"
                value={dataLimit}
                onChange={(e) => setDataLimit(Math.max(1, parseInt(e.target.value) || 100000))}
                style={{
                  width: "100%",
                  padding: "4px 8px",
                  fontSize: 12,
                  border: "1px solid var(--border-color)",
                  borderRadius: 4,
                }}
              />
            </div>

            <button
              onClick={handleDataExport}
              disabled={!state.info}
              style={{
                width: "100%",
                padding: "8px 16px",
                background: "var(--accent)",
                color: "white",
                border: "none",
                borderRadius: 4,
                fontSize: 12,
                fontWeight: 500,
                cursor: state.info ? "pointer" : "not-allowed",
                opacity: state.info ? 1 : 0.5,
              }}
            >
              Download {dataFormat.toUpperCase()}
            </button>
          </div>

          {/* Graph Export Section */}
          <div>
            <div
              style={{
                fontSize: 11,
                fontWeight: 600,
                color: "var(--text-secondary)",
                marginBottom: 8,
                textTransform: "uppercase",
                letterSpacing: 0.5,
              }}
            >
              Network Graph
            </div>

            {/* Format selector */}
            <div style={{ display: "flex", gap: 4, marginBottom: 8 }}>
              {graphFormats.map((fmt) => (
                <button
                  key={fmt.value}
                  onClick={() => setGraphFormat(fmt.value)}
                  title={fmt.desc}
                  style={{
                    flex: 1,
                    padding: "6px 4px",
                    fontSize: 10,
                    border: `1px solid ${graphFormat === fmt.value ? "var(--accent)" : "var(--border-color)"}`,
                    borderRadius: 4,
                    background: graphFormat === fmt.value ? "var(--accent)" : "white",
                    color: graphFormat === fmt.value ? "white" : "var(--text-primary)",
                    cursor: "pointer",
                    whiteSpace: "nowrap",
                    overflow: "hidden",
                  }}
                >
                  {fmt.label}
                </button>
              ))}
            </div>

            {/* Options */}
            <div style={{ display: "flex", flexDirection: "column", gap: 6, marginBottom: 8 }}>
              <label style={{ display: "flex", alignItems: "center", gap: 8, cursor: "pointer" }}>
                <input
                  type="checkbox"
                  checked={includeCommunities}
                  onChange={(e) => setIncludeCommunities(e.target.checked)}
                  style={{ width: 14, height: 14 }}
                />
                <span style={{ fontSize: 11 }}>Include community assignments</span>
              </label>
              <label style={{ display: "flex", alignItems: "center", gap: 8, cursor: "pointer" }}>
                <input
                  type="checkbox"
                  checked={includeLayout}
                  onChange={(e) => setIncludeLayout(e.target.checked)}
                  style={{ width: 14, height: 14 }}
                />
                <span style={{ fontSize: 11 }}>Include node positions</span>
              </label>
            </div>

            <button
              onClick={handleGraphExport}
              disabled={!state.info}
              style={{
                width: "100%",
                padding: "8px 16px",
                background: "var(--accent)",
                color: "white",
                border: "none",
                borderRadius: 4,
                fontSize: 12,
                fontWeight: 500,
                cursor: state.info ? "pointer" : "not-allowed",
                opacity: state.info ? 1 : 0.5,
              }}
            >
              Download Graph
            </button>
          </div>

          {/* Current filter info */}
          <div
            style={{
              padding: 8,
              background: "var(--bg-primary)",
              borderRadius: 4,
              fontSize: 10,
              color: "var(--text-secondary)",
            }}
          >
            Export will use current filter: <strong>{state.metric}</strong> &ge; {state.cutoff}
          </div>
        </div>
      )}
    </div>
  );
}
