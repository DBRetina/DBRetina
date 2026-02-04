import { useEffect, useCallback } from "react";
import { useDashboard } from "./state/context";
import { fetchInfo, fetchGraphData, fetchNeighborhood, getAPIError, isAPIError } from "./api";
import Sidebar from "./components/Sidebar/Sidebar";
import Toolbar from "./components/MainPanel/Toolbar";
import NetworkView from "./components/MainPanel/NetworkView";
import TableView from "./components/MainPanel/TableView";
import StatsView from "./components/MainPanel/StatsView";
import QueryPanel from "./components/QueryPanel/QueryPanel";
import DetailPanel from "./components/DetailPanel/DetailPanel";
import { ErrorBoundary, LoadingState } from "./components/Common";

export default function App() {
  const { state, dispatch } = useDashboard();

  // Handle errors with structured API error support
  const handleError = useCallback(
    (error: unknown, category: "global" | "graph" | "query" = "global") => {
      const apiError = getAPIError(error);
      if (apiError) {
        dispatch({ type: "SET_API_ERROR", category, error: apiError });
      } else {
        // Fallback for non-API errors
        const message = error instanceof Error ? error.message : String(error);
        dispatch({
          type: "SET_API_ERROR",
          category,
          error: {
            error: true,
            error_code: "CLIENT_ERROR",
            detail: message,
          },
        });
      }
    },
    [dispatch]
  );

  // Clear specific error
  const clearError = useCallback(
    (category: "global" | "graph" | "query") => {
      dispatch({ type: "SET_API_ERROR", category, error: null });
    },
    [dispatch]
  );

  // Load dataset info on mount
  useEffect(() => {
    dispatch({ type: "START_LOADING", key: "info", message: "Loading dataset..." });

    fetchInfo()
      .then((info) => {
        dispatch({ type: "SET_INFO", info });
        clearError("global");
      })
      .catch((e) => handleError(e, "global"))
      .finally(() => dispatch({ type: "STOP_LOADING", key: "info" }));
  }, [dispatch, handleError, clearError]);

  // Fetch graph data when metric/cutoff/focus changes
  useEffect(() => {
    if (!state.info) return;

    dispatch({ type: "START_LOADING", key: "graph", message: "Loading graph..." });

    const promise = state.focusGroup
      ? fetchNeighborhood(state.focusGroup, state.metric, state.cutoff)
      : fetchGraphData(state.metric, state.cutoff);

    promise
      .then((data) => {
        dispatch({ type: "SET_GRAPH_DATA", data });
        clearError("graph");

        // Track truncation state
        dispatch({
          type: "SET_DATA_STATE",
          state: {
            isEmpty: data.nodes.length === 0,
            isTruncated: data.meta.truncated,
            truncationInfo: data.meta.truncated
              ? { returned: data.meta.returned_nodes, total: data.meta.total_nodes }
              : undefined,
          },
        });
      })
      .catch((e) => handleError(e, "graph"))
      .finally(() => dispatch({ type: "STOP_LOADING", key: "graph" }));
  }, [state.info, state.metric, state.cutoff, state.focusGroup, dispatch, handleError, clearError]);

  // Get current loading state
  const graphLoading = state.loadingProgress.graph;
  const infoLoading = state.loadingProgress.info;

  // Get current errors
  const globalError = state.errors.global;
  const graphError = state.errors.graph;

  // Render error bar with dismiss button
  const renderErrorBar = () => {
    const error = globalError || graphError;
    if (!error) return null;

    return (
      <div className="error-bar" style={{ display: "flex", alignItems: "center", gap: 8 }}>
        <span style={{ flex: 1 }}>
          {error.detail}
          {error.context?.suggestion ? (
            <span style={{ opacity: 0.8, marginLeft: 8 }}>
              Suggestion: {String(error.context.suggestion)}
            </span>
          ) : null}
        </span>
        <button
          onClick={() => {
            if (globalError) clearError("global");
            if (graphError) clearError("graph");
          }}
          style={{
            border: "none",
            background: "transparent",
            color: "inherit",
            cursor: "pointer",
            padding: "2px 6px",
            fontSize: 16,
          }}
        >
          &times;
        </button>
      </div>
    );
  };

  // Render truncation warning
  const renderTruncationWarning = () => {
    if (!state.dataState.isTruncated || !state.dataState.truncationInfo) return null;

    const { returned, total } = state.dataState.truncationInfo;
    return (
      <div
        style={{
          padding: "6px 16px",
          background: "var(--warning)",
          color: "#000",
          fontSize: 12,
          display: "flex",
          alignItems: "center",
          gap: 8,
        }}
      >
        <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
          <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm1 15h-2v-2h2v2zm0-4h-2V7h2v6z" />
        </svg>
        <span>
          Showing {returned.toLocaleString()} of {total.toLocaleString()} nodes.
          Increase cutoff to reduce graph size.
        </span>
      </div>
    );
  };

  return (
    <ErrorBoundary
      fallback={
        <div style={{ padding: 48, textAlign: "center" }}>
          <h2>Something went wrong</h2>
          <p>The dashboard encountered an unexpected error.</p>
          <button
            onClick={() => window.location.reload()}
            style={{
              marginTop: 16,
              padding: "8px 16px",
              background: "var(--accent)",
              color: "white",
              border: "none",
              borderRadius: 4,
              cursor: "pointer",
            }}
          >
            Reload Dashboard
          </button>
        </div>
      }
    >
      <div className="app-shell">
        <Sidebar />
        <div className="main-panel">
          {renderErrorBar()}
          {renderTruncationWarning()}
          <Toolbar />
          <div className="content-area">
            {/* Enhanced loading overlay with message */}
            {(graphLoading?.active || infoLoading?.active) && (
              <LoadingState
                message={graphLoading?.message || infoLoading?.message}
                overlay
              />
            )}
            {state.activeView === "network" && <NetworkView />}
            {state.activeView === "table" && <TableView />}
            {state.activeView === "stats" && <StatsView />}
            {(state.selectedNode || state.selectedEdge) && <DetailPanel />}
          </div>
          {state.queryPanelOpen && <QueryPanel />}
        </div>
      </div>
    </ErrorBoundary>
  );
}
