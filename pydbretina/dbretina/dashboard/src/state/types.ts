export interface DatasetInfo {
  num_pairs: number;
  num_groups: number;
  has_pvalue: boolean;
  has_genes: boolean;
  default_metric: string;
  default_cutoff: number;
  valid_metrics: string[];
  manifest: Record<string, unknown>;
  statistics?: {
    [metric: string]: {
      histogram: number[];
      min: number;
      max: number;
      mean: number;
    };
  };
}

export interface GraphNode {
  id: string;
  label: string;
  degree: number;
  community: number;
  pagerank: number;
}

export interface GraphEdge {
  source: string;
  target: string;
  weight: number;
  shared_features: number;
}

export interface GraphData {
  nodes: GraphNode[];
  edges: GraphEdge[];
  meta: {
    total_nodes: number;
    total_edges: number;
    returned_nodes: number;
    returned_edges: number;
    truncated: boolean;
    metric: string;
    cutoff: number;
  };
}

export interface QueryResult {
  columns: string[];
  row_count: number;
  rows: Record<string, unknown>[];
}

export interface GroupMatch {
  id: number;
  name: string;
}

/** Structured API error from backend */
export interface APIError {
  error: boolean;
  error_code: string;
  detail: string;
  context?: Record<string, unknown>;
}

/** Loading state with optional progress and timing */
export interface LoadingProgress {
  active: boolean;
  message: string;
  progress?: number; // 0-100 for determinate progress
  startTime?: number; // timestamp for timeout warnings
}

/** Compare mode state for two-node analysis */
export interface CompareState {
  active: boolean;
  nodeA: GraphNode | null;
  nodeB: GraphNode | null;
}

/** Metric profile for a single metric */
export interface MetricProfile {
  metric: string;
  avg: number;
  max: number;
  count: number;
}

/** Data state indicators for UI feedback */
export interface DataState {
  isEmpty: boolean;
  isTruncated: boolean;
  truncationInfo?: {
    returned: number;
    total: number;
  };
}

export type LayoutAlgorithm = "force" | "fr" | "fa2" | "drl" | "kk" | "circle" | "grid";

/** Graph renderer type */
export type RendererType = "forcegraph" | "cytoscape" | "sigma" | "auto";

/** Selection mode for graph interactions */
export type SelectionMode = "click" | "lasso" | "box";

/** Selection state for multi-node/edge selection */
export interface SelectionState {
  mode: SelectionMode;
  selectedNodes: Set<string>;
  selectedEdges: Set<string>;
}

/** Gene importance score from backend */
export interface GeneScore {
  gene: string;
  score: number;
  global_freq?: number;
  specificity?: number;
  num_diseases?: number;
}

/** Result from hub genes analysis */
export interface HubGenesResult {
  group: string;
  method: string;
  hops: number;
  genes: GeneScore[];
}

/** Result from explain pair analysis */
export interface ExplainPairResult {
  group_a: string;
  group_b: string;
  method: string;
  gene_count: number;
  genes: GeneScore[];
}

/** Shared genes along a path edge */
export interface PathEdgeGenes {
  from: string;
  to: string;
  shared_count: number;
  genes: string[];
}

/** Result from shortest path query */
export interface PathResult {
  source: string;
  target: string;
  path_length: number;
  path_nodes: string[];
  connected: boolean;
  shared_genes_along_path?: PathEdgeGenes[];
}

/** Gene statistics for a single gene */
export interface GeneStatistics {
  gene: string;
  group_count: number;
  total_groups: number;
  prevalence_percent: number;
}

/** Group info for gene groups query */
export interface GeneGroupInfo {
  id: number;
  name: string;
}

/** Result from gene groups query */
export interface GeneGroupsResult {
  gene: string;
  group_count: number;
  groups: GeneGroupInfo[];
}

export interface DashboardState {
  info: DatasetInfo | null;
  metric: string;
  cutoff: number;
  graphData: GraphData | null;
  activeView: "network" | "table" | "stats";
  selectedNode: GraphNode | null;
  selectedEdge: { source: string; target: string } | null;
  queryPanelOpen: boolean;
  focusGroup: string | null;
  highlightNodes: Set<string>;

  // Layout settings
  layoutAlgorithm: LayoutAlgorithm;
  layoutPositions: Record<string, [number, number]> | null;

  // Renderer settings
  rendererType: RendererType;

  // Selection state
  selection: SelectionState;

  // Enhanced error tracking
  errors: {
    global: APIError | null;
    graph: APIError | null;
    query: APIError | null;
  };

  // Enhanced loading with progress
  loadingProgress: Record<string, LoadingProgress>;

  // Data state tracking
  dataState: DataState;

  // Biological analysis state
  hubGenesResult: HubGenesResult | null;
  explainPairResult: ExplainPairResult | null;
  pathResult: PathResult | null;
  activeGene: string | null;
  geneStatistics: GeneStatistics | null;
  geneGroups: GeneGroupsResult | null;

  // Compare mode
  compareState: CompareState;

  // Detail panel visibility
  detailPanelOpen: boolean;

  // Node filter — when set, only these nodes (and edges between them) are shown
  nodeFilter: Set<string> | null;

  // Legacy fields (kept for backwards compatibility)
  loading: Record<string, boolean>;
  error: string | null;
}

export type DashboardAction =
  | { type: "SET_INFO"; info: DatasetInfo }
  | { type: "SET_METRIC"; metric: string }
  | { type: "SET_CUTOFF"; cutoff: number }
  | { type: "SET_GRAPH_DATA"; data: GraphData }
  | { type: "SET_ACTIVE_VIEW"; view: "network" | "table" | "stats" }
  | { type: "SELECT_NODE"; node: GraphNode | null }
  | { type: "SELECT_EDGE"; edge: { source: string; target: string } | null }
  | { type: "TOGGLE_QUERY_PANEL" }
  | { type: "SET_FOCUS_GROUP"; group: string | null }
  | { type: "SET_HIGHLIGHT_NODES"; nodes: Set<string> }
  | { type: "SET_LOADING"; key: string; value: boolean }
  | { type: "SET_ERROR"; error: string | null }
  // Enhanced error actions
  | { type: "SET_API_ERROR"; category: "global" | "graph" | "query"; error: APIError | null }
  | { type: "CLEAR_ALL_ERRORS" }
  // Enhanced loading actions
  | { type: "START_LOADING"; key: string; message?: string }
  | { type: "UPDATE_LOADING_PROGRESS"; key: string; progress: number }
  | { type: "STOP_LOADING"; key: string }
  // Data state actions
  | { type: "SET_DATA_STATE"; state: Partial<DataState> }
  // Layout actions
  | { type: "SET_LAYOUT_ALGORITHM"; algorithm: LayoutAlgorithm }
  | { type: "SET_LAYOUT_POSITIONS"; positions: Record<string, [number, number]> | null }
  // Renderer actions
  | { type: "SET_RENDERER_TYPE"; renderer: RendererType }
  // Selection actions
  | { type: "SET_SELECTION_MODE"; mode: SelectionMode }
  | { type: "ADD_SELECTED_NODES"; nodes: string[] }
  | { type: "REMOVE_SELECTED_NODES"; nodes: string[] }
  | { type: "SET_SELECTED_NODES"; nodes: Set<string> }
  | { type: "CLEAR_SELECTION" }
  | { type: "INVERT_SELECTION" }
  // Biological analysis actions
  | { type: "SET_HUB_GENES"; result: HubGenesResult | null }
  | { type: "SET_EXPLAIN_PAIR"; result: ExplainPairResult | null }
  | { type: "SET_PATH_RESULT"; result: PathResult | null }
  | { type: "SET_ACTIVE_GENE"; gene: string | null }
  | { type: "SET_GENE_STATISTICS"; stats: GeneStatistics | null }
  | { type: "SET_GENE_GROUPS"; result: GeneGroupsResult | null }
  | { type: "CLEAR_BIOLOGICAL_ANALYSIS" }
  // Compare mode actions
  | { type: "ENTER_COMPARE_MODE"; nodeA: GraphNode }
  | { type: "SET_COMPARE_NODE_B"; nodeB: GraphNode }
  | { type: "EXIT_COMPARE_MODE" }
  // Detail panel actions
  | { type: "TOGGLE_DETAIL_PANEL" }
  | { type: "SET_DETAIL_PANEL_OPEN"; open: boolean }
  // Node filter actions
  | { type: "SET_NODE_FILTER"; nodes: Set<string> }
  | { type: "CLEAR_NODE_FILTER" };
