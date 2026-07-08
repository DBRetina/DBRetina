import { DashboardState, DashboardAction, LoadingProgress } from "./types";

export const initialState: DashboardState = {
  info: null,
  metric: "ochiai",
  cutoff: 20,
  graphData: null,
  activeView: "network",
  selectedNode: null,
  selectedEdge: null,
  queryPanelOpen: false,
  focusGroup: null,
  highlightNodes: new Set(),

  // Layout settings
  layoutAlgorithm: "force",
  layoutPositions: null,

  // Renderer settings
  rendererType: "forcegraph",

  // Selection state
  selection: {
    mode: "click",
    selectedNodes: new Set(),
    selectedEdges: new Set(),
  },

  // Enhanced error tracking
  errors: {
    global: null,
    graph: null,
    query: null,
  },

  // Enhanced loading with progress
  loadingProgress: {},

  // Data state tracking
  dataState: {
    isEmpty: false,
    isTruncated: false,
  },

  // Biological analysis state
  hubGenesResult: null,
  explainPairResult: null,
  pathResult: null,
  activeGene: null,
  geneStatistics: null,
  geneGroups: null,

  // Compare mode
  compareState: {
    active: false,
    nodeA: null,
    nodeB: null,
  },

  // Detail panel starts collapsed
  detailPanelOpen: false,

  // Node filter
  nodeFilter: null,

  // Edge filter (Advanced Filter)
  edgeFilter: null,

  // Legacy fields
  loading: {},
  error: null,
};

export function dashboardReducer(
  state: DashboardState,
  action: DashboardAction
): DashboardState {
  switch (action.type) {
    case "SET_INFO":
      return {
        ...state,
        info: action.info,
        metric: action.info.default_metric,
        cutoff: action.info.default_cutoff,
      };
    case "SET_METRIC":
      return { ...state, metric: action.metric, graphData: null };
    case "SET_CUTOFF":
      return { ...state, cutoff: action.cutoff, graphData: null };
    case "SET_GRAPH_DATA":
      return {
        ...state,
        graphData: action.data,
        // Exit compare mode and clear filters when graph data changes
        compareState: { active: false, nodeA: null, nodeB: null },
        nodeFilter: null,
        edgeFilter: null,
      };
    case "SET_ACTIVE_VIEW":
      return { ...state, activeView: action.view };
    case "SELECT_NODE":
      return {
        ...state,
        selectedNode: action.node,
        selectedEdge: null,
        // Auto-open detail panel when selecting a node
        detailPanelOpen: action.node ? true : state.detailPanelOpen,
      };
    case "SELECT_EDGE":
      return {
        ...state,
        selectedEdge: action.edge,
        selectedNode: null,
        // Auto-open detail panel when selecting an edge
        detailPanelOpen: action.edge ? true : state.detailPanelOpen,
      };
    case "TOGGLE_QUERY_PANEL":
      return { ...state, queryPanelOpen: !state.queryPanelOpen };
    case "SET_FOCUS_GROUP":
      return { ...state, focusGroup: action.group };
    case "SET_HIGHLIGHT_NODES":
      return { ...state, highlightNodes: action.nodes };
    case "SET_LOADING":
      return {
        ...state,
        loading: { ...state.loading, [action.key]: action.value },
      };
    case "SET_ERROR":
      return { ...state, error: action.error };

    // Enhanced error handling
    case "SET_API_ERROR":
      return {
        ...state,
        errors: { ...state.errors, [action.category]: action.error },
        // Also set legacy error for backwards compatibility
        error: action.error?.detail || null,
      };
    case "CLEAR_ALL_ERRORS":
      return {
        ...state,
        errors: { global: null, graph: null, query: null },
        error: null,
      };

    // Enhanced loading with progress
    case "START_LOADING": {
      const loadingState: LoadingProgress = {
        active: true,
        message: action.message || "Loading...",
        startTime: Date.now(),
      };
      return {
        ...state,
        loadingProgress: { ...state.loadingProgress, [action.key]: loadingState },
        loading: { ...state.loading, [action.key]: true },
      };
    }
    case "UPDATE_LOADING_PROGRESS":
      if (state.loadingProgress[action.key]) {
        return {
          ...state,
          loadingProgress: {
            ...state.loadingProgress,
            [action.key]: {
              ...state.loadingProgress[action.key],
              progress: action.progress,
            },
          },
        };
      }
      return state;
    case "STOP_LOADING": {
      const { [action.key]: _, ...remainingProgress } = state.loadingProgress;
      return {
        ...state,
        loadingProgress: remainingProgress,
        loading: { ...state.loading, [action.key]: false },
      };
    }

    // Data state
    case "SET_DATA_STATE":
      return {
        ...state,
        dataState: { ...state.dataState, ...action.state },
      };

    // Layout actions
    case "SET_LAYOUT_ALGORITHM":
      return {
        ...state,
        layoutAlgorithm: action.algorithm,
        layoutPositions: null, // Clear positions when algorithm changes
      };
    case "SET_LAYOUT_POSITIONS":
      return {
        ...state,
        layoutPositions: action.positions,
      };

    // Renderer actions
    case "SET_RENDERER_TYPE":
      return {
        ...state,
        rendererType: action.renderer,
      };

    // Selection actions
    case "SET_SELECTION_MODE":
      return {
        ...state,
        selection: { ...state.selection, mode: action.mode },
      };
    case "ADD_SELECTED_NODES": {
      const newSelected = new Set(state.selection.selectedNodes);
      for (const nodeId of action.nodes) {
        newSelected.add(nodeId);
      }
      return {
        ...state,
        selection: { ...state.selection, selectedNodes: newSelected },
      };
    }
    case "REMOVE_SELECTED_NODES": {
      const newSelected = new Set(state.selection.selectedNodes);
      for (const nodeId of action.nodes) {
        newSelected.delete(nodeId);
      }
      return {
        ...state,
        selection: { ...state.selection, selectedNodes: newSelected },
      };
    }
    case "SET_SELECTED_NODES":
      return {
        ...state,
        selection: { ...state.selection, selectedNodes: action.nodes },
      };
    case "CLEAR_SELECTION":
      return {
        ...state,
        selection: {
          ...state.selection,
          selectedNodes: new Set(),
          selectedEdges: new Set(),
        },
      };
    case "INVERT_SELECTION": {
      if (!state.graphData) return state;
      const allNodeIds = new Set(state.graphData.nodes.map((n) => n.id));
      const inverted = new Set<string>();
      for (const nodeId of allNodeIds) {
        if (!state.selection.selectedNodes.has(nodeId)) {
          inverted.add(nodeId);
        }
      }
      return {
        ...state,
        selection: { ...state.selection, selectedNodes: inverted },
      };
    }

    // Biological analysis actions
    case "SET_HUB_GENES":
      return { ...state, hubGenesResult: action.result };
    case "SET_EXPLAIN_PAIR":
      return { ...state, explainPairResult: action.result };
    case "SET_PATH_RESULT":
      return { ...state, pathResult: action.result };
    case "SET_ACTIVE_GENE":
      return { ...state, activeGene: action.gene };
    case "SET_GENE_STATISTICS":
      return { ...state, geneStatistics: action.stats };
    case "SET_GENE_GROUPS":
      return { ...state, geneGroups: action.result };
    case "CLEAR_BIOLOGICAL_ANALYSIS":
      return {
        ...state,
        hubGenesResult: null,
        explainPairResult: null,
        pathResult: null,
        activeGene: null,
        geneStatistics: null,
        geneGroups: null,
      };

    // Compare mode actions
    case "ENTER_COMPARE_MODE":
      return {
        ...state,
        compareState: { active: true, nodeA: action.nodeA, nodeB: null },
        selectedNode: null,
        selectedEdge: null,
        detailPanelOpen: true,
      };
    case "SET_COMPARE_NODE_B":
      return {
        ...state,
        compareState: { ...state.compareState, nodeB: action.nodeB },
      };
    case "EXIT_COMPARE_MODE":
      return {
        ...state,
        compareState: { active: false, nodeA: null, nodeB: null },
      };

    // Detail panel actions
    case "TOGGLE_DETAIL_PANEL":
      return { ...state, detailPanelOpen: !state.detailPanelOpen };
    case "SET_DETAIL_PANEL_OPEN":
      return { ...state, detailPanelOpen: action.open };

    // Node filter actions
    case "SET_NODE_FILTER":
      return { ...state, nodeFilter: action.nodes };
    case "CLEAR_NODE_FILTER":
      return { ...state, nodeFilter: null };

    // Edge filter actions (Advanced Filter)
    case "SET_EDGE_FILTER":
      return { ...state, edgeFilter: action.edges };
    case "CLEAR_EDGE_FILTER":
      return { ...state, edgeFilter: null };

    default:
      return state;
  }
}
