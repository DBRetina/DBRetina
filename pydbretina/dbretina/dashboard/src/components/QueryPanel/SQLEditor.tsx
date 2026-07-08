import { useState, useCallback, useMemo } from "react";
import CodeMirror from "@uiw/react-codemirror";
import { sql } from "@codemirror/lang-sql";
import { oneDark } from "@codemirror/theme-one-dark";
import { keymap } from "@codemirror/view";
import { executeSql, getAPIError } from "../../api";
import { QueryResult, APIError } from "../../state/types";
import { validateSQLSafety } from "../../utils/validation";
import { useQueryHistory, SQL_TEMPLATES } from "../../hooks/useQueryHistory";

interface Props {
  onResult: (result: QueryResult) => void;
  onError: (error: string | APIError) => void;
  onLoading?: (loading: boolean) => void;
}

const SQL_PLACEHOLDER = `-- Example: Find top similar pairs
SELECT group_1_name, group_2_name, ochiai
FROM pairs
WHERE ochiai > 80
ORDER BY ochiai DESC
LIMIT 20`;

export default function SQLEditor({ onResult, onError, onLoading }: Props) {
  const [query, setQuery] = useState(SQL_PLACEHOLDER);
  const [running, setRunning] = useState(false);
  const [showHistory, setShowHistory] = useState(false);
  const [showTemplates, setShowTemplates] = useState(false);

  const { history, addToHistory } = useQueryHistory("sql");

  const handleRun = useCallback(async () => {
    const q = query.trim();
    if (!q) return;

    // Frontend validation
    const validation = validateSQLSafety(q);
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
      const result = await executeSql(q);
      onResult(result);
      addToHistory(q, true);
    } catch (e: unknown) {
      const apiError = getAPIError(e);
      if (apiError) {
        onError(apiError);
      } else {
        onError(e instanceof Error ? e.message : "SQL execution failed");
      }
      addToHistory(q, false);
    } finally {
      setRunning(false);
      onLoading?.(false);
    }
  }, [query, onResult, onError, onLoading, addToHistory]);

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
              {SQL_TEMPLATES.map((tpl) => (
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
        extensions={[sql(), runKeymap]}
        onChange={(value) => setQuery(value)}
        basicSetup={{ lineNumbers: false, foldGutter: false }}
      />
      <button className="query-run-btn" onClick={handleRun} disabled={running}>
        {running ? "Running..." : "Run (Ctrl+Enter)"}
      </button>
    </div>
  );
}
