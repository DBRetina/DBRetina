import { useState, useCallback } from "react";
import { useDashboard } from "../../state/context";
import SQLEditor from "./SQLEditor";
import CypherEditor from "./CypherEditor";
import ResultsTable from "./ResultsTable";
import { QueryResult, APIError } from "../../state/types";
import { EmptyState } from "../Common";

/**
 * Extract all string values from query results and find matching graph nodes.
 * Returns node IDs that match any string value in the results.
 */
function extractMatchingNodeIds(
  result: QueryResult,
  graphNodes: Array<{ id: string; label: string }>
): Set<string> {
  // Build a lookup: lowercase label → node id
  const labelToId = new Map<string, string>();
  for (const node of graphNodes) {
    labelToId.set(node.label.toLowerCase(), node.id);
  }

  const matched = new Set<string>();
  for (const row of result.rows) {
    for (const val of Object.values(row)) {
      if (typeof val === "string") {
        const nodeId = labelToId.get(val.toLowerCase());
        if (nodeId) {
          matched.add(nodeId);
        }
      }
    }
  }
  return matched;
}

export default function QueryPanel() {
  const { state, dispatch } = useDashboard();
  const [activeTab, setActiveTab] = useState<"sql" | "cypher">("sql");
  const [results, setResults] = useState<QueryResult | null>(null);
  const [error, setError] = useState<APIError | null>(null);
  const [loading, setLoading] = useState(false);
  const [highlightedCount, setHighlightedCount] = useState(0);

  const isOpen = state.queryPanelOpen;

  const handleResult = useCallback(
    (result: QueryResult) => {
      setResults(result);
      setError(null);

      // Auto-highlight matching nodes in the graph
      if (state.graphData && result.row_count > 0) {
        const matchedIds = extractMatchingNodeIds(result, state.graphData.nodes);
        if (matchedIds.size > 0) {
          dispatch({ type: "SET_HIGHLIGHT_NODES", nodes: matchedIds });
          setHighlightedCount(matchedIds.size);
        } else {
          setHighlightedCount(0);
        }
      } else {
        setHighlightedCount(0);
      }
    },
    [state.graphData, dispatch]
  );

  function handleError(err: string | APIError) {
    if (typeof err === "string") {
      setError({ error: true, error_code: "QUERY_ERROR", detail: err });
    } else {
      setError(err);
    }
    setResults(null);
    setHighlightedCount(0);
  }

  function handleLoading(isLoading: boolean) {
    setLoading(isLoading);
  }

  function handleClearHighlights() {
    dispatch({ type: "SET_HIGHLIGHT_NODES", nodes: new Set() });
    setHighlightedCount(0);
  }

  function handleFocusNodes() {
    if (state.highlightNodes.size > 0) {
      dispatch({ type: "SET_NODE_FILTER", nodes: new Set(state.highlightNodes) });
    }
  }

  function handleResetFilter() {
    dispatch({ type: "CLEAR_NODE_FILTER" });
  }

  function handleToggle() {
    dispatch({ type: "TOGGLE_QUERY_PANEL" });
  }

  // Determine error styling based on error type
  const getErrorStyle = () => {
    if (!error) return {};
    if (error.error_code === "UNSAFE_QUERY") {
      return { background: "#fff3cd", borderColor: "#ffc107" };
    }
    return {};
  };

  const hasResults = !!(results && results.row_count > 0);

  return (
    <div className={`query-panel-wrapper ${isOpen ? "open" : "collapsed"}`}>
      {/* Toggle tab — always visible */}
      <button
        className="query-panel-toggle"
        onClick={handleToggle}
        title={isOpen ? "Collapse query panel" : "Expand query panel"}
      >
        <svg width="8" height="14" viewBox="0 0 8 14" fill="none">
          {isOpen ? (
            <path d="M7 1L1 7l6 6" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
          ) : (
            <path d="M1 1l6 6-6 6" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
          )}
        </svg>
        {!isOpen && hasResults && <span className="query-panel-badge" />}
      </button>

      {/* Panel content */}
      <div className="query-panel">
        <div className="query-panel-header">
          <span>Query</span>
          <div className="query-tabs">
            <button
              className={`query-tab ${activeTab === "sql" ? "active" : ""}`}
              onClick={() => setActiveTab("sql")}
            >
              SQL
            </button>
            <button
              className={`query-tab ${activeTab === "cypher" ? "active" : ""}`}
              onClick={() => setActiveTab("cypher")}
            >
              Cypher
            </button>
          </div>
          {loading && (
            <span style={{ fontSize: 11, color: "var(--accent)" }}>
              Running...
            </span>
          )}
          {!loading && results && (
            <span style={{ fontSize: 11, color: "var(--text-secondary)" }}>
              {results.row_count} rows
            </span>
          )}
        </div>

        <div className="query-panel-body">
          {activeTab === "sql" && (
            <SQLEditor onResult={handleResult} onError={handleError} onLoading={handleLoading} />
          )}
          {activeTab === "cypher" && (
            <CypherEditor onResult={handleResult} onError={handleError} onLoading={handleLoading} />
          )}

          {/* Highlight indicator */}
          {highlightedCount > 0 && (
            <div
              style={{
                padding: "6px 12px",
                fontSize: 11,
                color: "var(--accent)",
                display: "flex",
                alignItems: "center",
                gap: 8,
                background: "#e8f0fe",
                margin: "0 12px",
                borderRadius: 4,
              }}
            >
              <svg width="12" height="12" viewBox="0 0 24 24" fill="currentColor">
                <circle cx="12" cy="12" r="10" />
              </svg>
              <span>{highlightedCount} highlighted</span>
              <div style={{ marginLeft: "auto", display: "flex", gap: 4 }}>
                {!state.nodeFilter ? (
                  <button
                    onClick={handleFocusNodes}
                    style={{
                      padding: "2px 8px",
                      fontSize: 10,
                      border: "1px solid var(--accent)",
                      borderRadius: 3,
                      background: "var(--accent)",
                      color: "white",
                      cursor: "pointer",
                    }}
                  >
                    Focus
                  </button>
                ) : (
                  <button
                    onClick={handleResetFilter}
                    style={{
                      padding: "2px 8px",
                      fontSize: 10,
                      border: "1px solid var(--danger)",
                      borderRadius: 3,
                      background: "var(--danger)",
                      color: "white",
                      cursor: "pointer",
                    }}
                  >
                    Reset
                  </button>
                )}
                <button
                  onClick={handleClearHighlights}
                  style={{
                    padding: "2px 8px",
                    fontSize: 10,
                    border: "1px solid var(--accent)",
                    borderRadius: 3,
                    background: "white",
                    color: "var(--accent)",
                    cursor: "pointer",
                  }}
                >
                  Clear
                </button>
              </div>
            </div>
          )}

          {/* Filter active indicator (shown when filter is set but no highlights) */}
          {state.nodeFilter && highlightedCount === 0 && (
            <div
              style={{
                padding: "6px 12px",
                fontSize: 11,
                color: "var(--warning, #e67e22)",
                display: "flex",
                alignItems: "center",
                gap: 8,
                background: "#fef3e2",
                margin: "0 12px",
                borderRadius: 4,
              }}
            >
              <svg width="12" height="12" viewBox="0 0 24 24" fill="currentColor">
                <path d="M10 18h4v-2h-4v2zM3 6v2h18V6H3zm3 7h12v-2H6v2z" />
              </svg>
              <span>{state.nodeFilter.size} filtered</span>
              <button
                onClick={handleResetFilter}
                style={{
                  marginLeft: "auto",
                  padding: "2px 8px",
                  fontSize: 10,
                  border: "1px solid var(--danger)",
                  borderRadius: 3,
                  background: "var(--danger)",
                  color: "white",
                  cursor: "pointer",
                }}
              >
                Reset
              </button>
            </div>
          )}

          {error && (
            <div
              style={{
                padding: "8px 12px",
                color: error.error_code === "UNSAFE_QUERY" ? "#856404" : "var(--danger)",
                fontSize: 12,
                display: "flex",
                alignItems: "flex-start",
                gap: 8,
                borderRadius: 4,
                margin: "0 12px",
                ...getErrorStyle(),
              }}
            >
              <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor" style={{ flexShrink: 0, marginTop: 1 }}>
                <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm1 15h-2v-2h2v2zm0-4h-2V7h2v6z" />
              </svg>
              <div>
                <div style={{ fontWeight: 500 }}>{error.detail}</div>
                {error.context?.suggestion ? (
                  <div style={{ fontSize: 11, opacity: 0.8, marginTop: 2 }}>
                    {String(error.context.suggestion)}
                  </div>
                ) : null}
                {error.error_code === "UNSAFE_QUERY" && (
                  <div style={{ fontSize: 11, opacity: 0.8, marginTop: 2 }}>
                    Only read queries (SELECT, MATCH, RETURN) are allowed.
                  </div>
                )}
              </div>
            </div>
          )}
          {results && results.row_count === 0 && (
            <div style={{ padding: "16px 12px" }}>
              <EmptyState
                type="no-results"
                title="No results"
                message="The query returned no rows."
                compact
              />
            </div>
          )}
          {results && results.row_count > 0 && <ResultsTable result={results} />}
        </div>
      </div>
    </div>
  );
}
