import { useState, useCallback } from "react";
import { useDashboard } from "../../state/context";
import { fetchGeneStatistics, fetchGeneGroups } from "../../api";

export default function GeneExplorer() {
  const { state, dispatch } = useDashboard();
  const [expanded, setExpanded] = useState(false);
  const [searchQuery, setSearchQuery] = useState("");
  const [searching, setSearching] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleSearch = useCallback(async () => {
    const query = searchQuery.trim();
    if (!query || query.length < 2) {
      dispatch({ type: "SET_GENE_STATISTICS", stats: null });
      dispatch({ type: "SET_GENE_GROUPS", result: null });
      return;
    }

    setSearching(true);
    setError(null);

    try {
      const [stats, groups] = await Promise.all([
        fetchGeneStatistics(query),
        fetchGeneGroups(query, 50),
      ]);

      dispatch({ type: "SET_GENE_STATISTICS", stats });
      dispatch({ type: "SET_GENE_GROUPS", result: groups });
      dispatch({ type: "SET_ACTIVE_GENE", gene: stats.gene });
    } catch (err) {
      setError(err instanceof Error ? err.message : "Gene not found");
      dispatch({ type: "SET_GENE_STATISTICS", stats: null });
      dispatch({ type: "SET_GENE_GROUPS", result: null });
    } finally {
      setSearching(false);
    }
  }, [searchQuery, dispatch]);

  const handleHighlightAll = useCallback(() => {
    if (!state.geneGroups || !state.graphData) return;

    const groupNames = new Set(state.geneGroups.groups.map((g) => g.name.toLowerCase()));
    const nodeIds = state.graphData.nodes
      .filter((n) => groupNames.has(n.label.toLowerCase()))
      .map((n) => n.id);

    dispatch({ type: "SET_HIGHLIGHT_NODES", nodes: new Set(nodeIds) });
  }, [state.geneGroups, state.graphData, dispatch]);

  const handleClearHighlight = useCallback(() => {
    dispatch({ type: "SET_HIGHLIGHT_NODES", nodes: new Set() });
  }, [dispatch]);

  const handleGroupClick = useCallback(
    (groupName: string) => {
      dispatch({ type: "SET_FOCUS_GROUP", group: groupName });
    },
    [dispatch]
  );

  const stats = state.geneStatistics;
  const groups = state.geneGroups;

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
        <label style={{ cursor: "pointer", marginBottom: 0 }}>Gene Explorer</label>
        <span style={{ fontSize: 12, color: "var(--text-secondary)" }}>
          {expanded ? "[-]" : "[+]"}
        </span>
      </div>

      {expanded && (
        <div style={{ display: "flex", flexDirection: "column", gap: 12 }}>
          {/* Search input - only searches on Enter */}
          <div style={{ display: "flex", gap: 4 }}>
            <input
              type="text"
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              onKeyDown={(e) => {
                if (e.key === "Enter") handleSearch();
              }}
              placeholder="Type gene name, press Enter..."
              style={{
                flex: 1,
                padding: "8px 12px",
                border: "1px solid var(--border-color)",
                borderRadius: 4,
                fontSize: 13,
              }}
            />
            <button
              onClick={handleSearch}
              disabled={searching}
              style={{
                padding: "8px 12px",
                background: "var(--accent)",
                color: "white",
                border: "none",
                borderRadius: 4,
                fontSize: 12,
                cursor: "pointer",
                whiteSpace: "nowrap",
              }}
            >
              {searching ? "..." : "Search"}
            </button>
          </div>
          {error && (
            <div style={{ fontSize: 11, color: "#e74c3c" }}>{error}</div>
          )}

          {/* Gene statistics */}
          {stats && (
            <div
              style={{
                padding: 12,
                background: "var(--bg-primary)",
                borderRadius: 6,
                border: "1px solid var(--border-color)",
              }}
            >
              <div
                style={{
                  fontSize: 14,
                  fontWeight: 600,
                  marginBottom: 8,
                  color: "var(--accent)",
                }}
              >
                {stats.gene}
              </div>
              <div style={{ display: "flex", gap: 16, fontSize: 12 }}>
                <div>
                  <div style={{ color: "var(--text-secondary)" }}>Groups</div>
                  <div style={{ fontWeight: 600, fontSize: 16 }}>{stats.group_count}</div>
                </div>
                <div>
                  <div style={{ color: "var(--text-secondary)" }}>Prevalence</div>
                  <div style={{ fontWeight: 600, fontSize: 16 }}>
                    {stats.prevalence_percent.toFixed(1)}%
                  </div>
                </div>
                <div>
                  <div style={{ color: "var(--text-secondary)" }}>Total</div>
                  <div style={{ fontWeight: 600, fontSize: 16 }}>{stats.total_groups}</div>
                </div>
              </div>
            </div>
          )}

          {/* Action buttons */}
          {groups && groups.groups.length > 0 && (
            <div style={{ display: "flex", gap: 8 }}>
              <button
                onClick={handleHighlightAll}
                style={{
                  flex: 1,
                  padding: "8px 12px",
                  background: "var(--accent)",
                  color: "white",
                  border: "none",
                  borderRadius: 4,
                  fontSize: 12,
                  cursor: "pointer",
                }}
              >
                Highlight All ({groups.groups.length})
              </button>
              {state.highlightNodes.size > 0 && (
                <button
                  onClick={handleClearHighlight}
                  style={{
                    padding: "8px 12px",
                    background: "white",
                    color: "var(--text-primary)",
                    border: "1px solid var(--border-color)",
                    borderRadius: 4,
                    fontSize: 12,
                    cursor: "pointer",
                  }}
                >
                  Clear
                </button>
              )}
            </div>
          )}

          {/* Groups list */}
          {groups && groups.groups.length > 0 && (
            <div>
              <div
                style={{
                  fontSize: 11,
                  fontWeight: 600,
                  color: "var(--text-secondary)",
                  marginBottom: 8,
                  textTransform: "uppercase",
                }}
              >
                Groups containing this gene
              </div>
              <div
                style={{
                  maxHeight: 200,
                  overflowY: "auto",
                  border: "1px solid var(--border-color)",
                  borderRadius: 4,
                }}
              >
                {groups.groups.map((group, i) => (
                  <div
                    key={group.id}
                    onClick={() => handleGroupClick(group.name)}
                    style={{
                      padding: "8px 12px",
                      fontSize: 12,
                      cursor: "pointer",
                      borderBottom:
                        i < groups.groups.length - 1 ? "1px solid var(--border-color)" : "none",
                      background: "white",
                    }}
                    onMouseOver={(e) => (e.currentTarget.style.background = "var(--bg-primary)")}
                    onMouseOut={(e) => (e.currentTarget.style.background = "white")}
                  >
                    {group.name}
                  </div>
                ))}
              </div>
            </div>
          )}

          {/* Help text */}
          {!stats && !error && (
            <div style={{ fontSize: 11, color: "var(--text-secondary)" }}>
              Enter a gene name (e.g., BRCA2, SRGAP2) and press Enter to see which
              groups contain it.
            </div>
          )}
        </div>
      )}
    </div>
  );
}
