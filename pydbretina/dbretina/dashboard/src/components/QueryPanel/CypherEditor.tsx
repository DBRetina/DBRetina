import { useState, useCallback, useMemo } from "react";
import CodeMirror from "@uiw/react-codemirror";
import { StreamLanguage } from "@codemirror/language";
import { oneDark } from "@codemirror/theme-one-dark";
import { keymap } from "@codemirror/view";
import { executeCypher, getAPIError } from "../../api";
import { useDashboard } from "../../state/context";
import { QueryResult, APIError } from "../../state/types";
import { validateCypherSafety } from "../../utils/validation";
import { useQueryHistory, CYPHER_TEMPLATES } from "../../hooks/useQueryHistory";

interface Props {
  onResult: (result: QueryResult) => void;
  onError: (error: string | APIError) => void;
  onLoading?: (loading: boolean) => void;
}

const CYPHER_PLACEHOLDER = `// Example: Find neighbors of a disease
MATCH (a:Group)-[r:SIMILAR_TO]->(b:Group)
WHERE a.name = "alzheimer's disease"
RETURN b.name, r.ochiai
ORDER BY r.ochiai DESC
LIMIT 20`;

// Simple Cypher syntax highlighting
const cypherLanguage = StreamLanguage.define({
  token(stream) {
    if (stream.match(/\/\/.*/)) return "comment";
    if (stream.match(/MATCH|WHERE|RETURN|ORDER BY|LIMIT|WITH|OPTIONAL|UNWIND|AS|AND|OR|NOT|IN|CONTAINS|STARTS WITH|ENDS WITH|IS NULL|IS NOT NULL|DESC|ASC/i)) {
      return "keyword";
    }
    if (stream.match(/:[A-Za-z_][A-Za-z0-9_]*/)) return "typeName";
    if (stream.match(/"[^"]*"|'[^']*'/)) return "string";
    if (stream.match(/\d+(\.\d+)?/)) return "number";
    if (stream.match(/[a-zA-Z_][a-zA-Z0-9_]*/)) return "variableName";
    stream.next();
    return null;
  },
});

export default function CypherEditor({ onResult, onError, onLoading }: Props) {
  const { state } = useDashboard();
  const [query, setQuery] = useState(CYPHER_PLACEHOLDER);
  const [running, setRunning] = useState(false);
  const [showHistory, setShowHistory] = useState(false);
  const [showTemplates, setShowTemplates] = useState(false);

  const { history, addToHistory } = useQueryHistory("cypher");

  const handleRun = useCallback(async () => {
    const q = query.trim();
    if (!q || q === CYPHER_PLACEHOLDER) return;

    // Frontend validation
    const validation = validateCypherSafety(q);
    if (!validation.valid) {
      onError({
        error: true,
        error_code: "UNSAFE_QUERY",
        detail: validation.error || "Invalid query",
        context: { blocked_operation: "frontend_validation" },
      });
      addToHistory(q, false);
      return;
    }

    setRunning(true);
    onLoading?.(true);

    try {
      const result = await executeCypher(q, state.metric, state.cutoff);
      onResult(result);
      addToHistory(q, true);
    } catch (e: unknown) {
      const apiError = getAPIError(e);
      if (apiError) {
        onError(apiError);
      } else {
        onError(e instanceof Error ? e.message : "Cypher execution failed");
      }
      addToHistory(q, false);
    } finally {
      setRunning(false);
      onLoading?.(false);
    }
  }, [query, state.metric, state.cutoff, onResult, onError, onLoading, addToHistory]);

  // Ctrl+Enter keymap extension
  const runKeymap = useMemo(
    () =>
      keymap.of([
        {
          key: "Ctrl-Enter",
          mac: "Cmd-Enter",
          run: () => {
            handleRun();
            return true;
          },
        },
      ]),
    [handleRun]
  );

  return (
    <div className="query-editor">
      {/* Toolbar */}
      <div style={{ display: "flex", gap: 4, marginBottom: 4 }}>
        {/* Templates dropdown */}
        <div style={{ position: "relative" }}>
          <button
            onClick={() => {
              setShowTemplates(!showTemplates);
              setShowHistory(false);
            }}
            style={{
              padding: "3px 8px",
              fontSize: 10,
              border: "1px solid var(--border-color)",
              borderRadius: 3,
              background: "white",
              cursor: "pointer",
            }}
          >
            Templates
          </button>
          {showTemplates && (
            <div
              style={{
                position: "absolute",
                top: "100%",
                left: 0,
                zIndex: 10,
                background: "white",
                border: "1px solid var(--border-color)",
                borderRadius: 4,
                boxShadow: "0 2px 8px rgba(0,0,0,0.15)",
                minWidth: 200,
                maxHeight: 200,
                overflow: "auto",
              }}
            >
              {CYPHER_TEMPLATES.map((tpl) => (
                <div
                  key={tpl.name}
                  onClick={() => {
                    setQuery(tpl.query);
                    setShowTemplates(false);
                  }}
                  style={{
                    padding: "6px 10px",
                    fontSize: 11,
                    cursor: "pointer",
                    borderBottom: "1px solid var(--border-color)",
                  }}
                  onMouseEnter={(e) => (e.currentTarget.style.background = "var(--bg-primary)")}
                  onMouseLeave={(e) => (e.currentTarget.style.background = "white")}
                >
                  {tpl.name}
                </div>
              ))}
            </div>
          )}
        </div>

        {/* History dropdown */}
        <div style={{ position: "relative" }}>
          <button
            onClick={() => {
              setShowHistory(!showHistory);
              setShowTemplates(false);
            }}
            style={{
              padding: "3px 8px",
              fontSize: 10,
              border: "1px solid var(--border-color)",
              borderRadius: 3,
              background: "white",
              cursor: "pointer",
            }}
          >
            History ({history.length})
          </button>
          {showHistory && history.length > 0 && (
            <div
              style={{
                position: "absolute",
                top: "100%",
                left: 0,
                zIndex: 10,
                background: "white",
                border: "1px solid var(--border-color)",
                borderRadius: 4,
                boxShadow: "0 2px 8px rgba(0,0,0,0.15)",
                minWidth: 300,
                maxHeight: 200,
                overflow: "auto",
              }}
            >
              {history.slice(0, 10).map((item, i) => (
                <div
                  key={i}
                  onClick={() => {
                    setQuery(item.query);
                    setShowHistory(false);
                  }}
                  style={{
                    padding: "6px 10px",
                    fontSize: 10,
                    fontFamily: "var(--font-mono)",
                    cursor: "pointer",
                    borderBottom: "1px solid var(--border-color)",
                    borderLeft: `3px solid ${item.success ? "var(--success)" : "var(--danger)"}`,
                    whiteSpace: "nowrap",
                    overflow: "hidden",
                    textOverflow: "ellipsis",
                  }}
                  onMouseEnter={(e) => (e.currentTarget.style.background = "var(--bg-primary)")}
                  onMouseLeave={(e) => (e.currentTarget.style.background = "white")}
                  title={item.query}
                >
                  {item.query.substring(0, 60)}...
                </div>
              ))}
            </div>
          )}
        </div>
      </div>

      <CodeMirror
        value={query}
        height="100px"
        theme={oneDark}
        extensions={[cypherLanguage, runKeymap]}
        onChange={(value) => setQuery(value)}
        basicSetup={{ lineNumbers: false, foldGutter: false }}
      />
      <button className="query-run-btn" onClick={handleRun} disabled={running}>
        {running ? "Running..." : "Run (Ctrl+Enter)"}
      </button>
    </div>
  );
}
