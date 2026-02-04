import { useState } from "react";
import SQLEditor from "./SQLEditor";
import CypherEditor from "./CypherEditor";
import ResultsTable from "./ResultsTable";
import { QueryResult, APIError } from "../../state/types";
import { EmptyState } from "../Common";

export default function QueryPanel() {
  const [activeTab, setActiveTab] = useState<"sql" | "cypher">("sql");
  const [results, setResults] = useState<QueryResult | null>(null);
  const [error, setError] = useState<APIError | null>(null);
  const [loading, setLoading] = useState(false);

  function handleResult(result: QueryResult) {
    setResults(result);
    setError(null);
  }

  function handleError(err: string | APIError) {
    if (typeof err === "string") {
      setError({ error: true, error_code: "QUERY_ERROR", detail: err });
    } else {
      setError(err);
    }
    setResults(null);
  }

  function handleLoading(isLoading: boolean) {
    setLoading(isLoading);
  }

  // Determine error styling based on error type
  const getErrorStyle = () => {
    if (!error) return {};
    if (error.error_code === "UNSAFE_QUERY") {
      return { background: "#fff3cd", borderColor: "#ffc107" };
    }
    return {};
  };

  return (
    <div className="query-panel">
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
        {loading && (
          <span style={{ marginLeft: "auto", fontSize: 11, color: "var(--accent)", padding: "6px 0" }}>
            Running...
          </span>
        )}
        {!loading && results && (
          <span style={{ marginLeft: "auto", fontSize: 11, color: "var(--text-secondary)", padding: "6px 0" }}>
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
  );
}
