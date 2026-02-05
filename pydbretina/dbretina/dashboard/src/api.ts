import { DatasetInfo, GraphData, QueryResult, GroupMatch, APIError } from "./state/types";

const BASE = "/api/v1";

/**
 * Custom error class that preserves structured API error information.
 */
export class APIRequestError extends Error {
  public readonly apiError: APIError | null;
  public readonly statusCode: number;

  constructor(message: string, statusCode: number, apiError: APIError | null = null) {
    super(message);
    this.name = "APIRequestError";
    this.statusCode = statusCode;
    this.apiError = apiError;
  }

  /** Check if this is a specific error code */
  isErrorCode(code: string): boolean {
    return this.apiError?.error_code === code;
  }

  /** Get structured error or fallback to generic */
  toAPIError(): APIError {
    return (
      this.apiError || {
        error: true,
        error_code: "UNKNOWN_ERROR",
        detail: this.message,
      }
    );
  }
}

async function handleResponse<T>(res: Response): Promise<T> {
  if (!res.ok) {
    let message = `API error ${res.status}`;
    let apiError: APIError | null = null;

    try {
      const body = await res.json();
      // Check if it's a structured error response
      if (body.error && body.error_code) {
        apiError = body as APIError;
        message = body.detail || message;
      } else {
        message = body.detail || body.error || body.message || message;
      }
    } catch {
      // If not JSON, try text but truncate
      try {
        const text = await res.text();
        if (text.length < 200 && !text.includes("<html")) {
          message = text;
        }
      } catch {
        // ignore
      }
    }

    throw new APIRequestError(message, res.status, apiError);
  }
  return res.json();
}

async function fetchJSON<T>(url: string, init?: RequestInit): Promise<T> {
  const res = await fetch(url, init);
  return handleResponse<T>(res);
}

export async function fetchInfo(): Promise<DatasetInfo> {
  return fetchJSON(`${BASE}/info`);
}

export async function fetchGraphData(
  metric: string,
  cutoff: number,
  limit = 5000
): Promise<GraphData> {
  return fetchJSON(
    `${BASE}/graph/data?metric=${metric}&cutoff=${cutoff}&limit=${limit}`
  );
}

export async function fetchNeighborhood(
  group: string,
  metric: string,
  cutoff: number,
  hops = 2
): Promise<GraphData> {
  return fetchJSON(
    `${BASE}/graph/neighborhood?group=${encodeURIComponent(group)}&metric=${metric}&cutoff=${cutoff}&hops=${hops}`
  );
}

export async function fetchLayout(
  metric: string,
  cutoff: number,
  algorithm = "fr"
): Promise<{ algorithm: string; positions: Record<string, [number, number]> }> {
  return fetchJSON(
    `${BASE}/graph/layout?metric=${metric}&cutoff=${cutoff}&algorithm=${algorithm}`
  );
}

export async function fetchCommunities(
  metric: string,
  cutoff: number,
  method = "leiden"
): Promise<{ communities: Record<string, number>; num_communities: number; sizes: Record<string, number> }> {
  return fetchJSON(
    `${BASE}/graph/communities?metric=${metric}&cutoff=${cutoff}&method=${method}`
  );
}

export async function searchGroups(q: string, limit = 20): Promise<{ matches: GroupMatch[] }> {
  return fetchJSON(`${BASE}/groups/search?q=${encodeURIComponent(q)}&limit=${limit}`);
}

export async function searchGenes(q: string, limit = 20): Promise<{ query: string; genes: Array<{ gene: string; group_count: number }> }> {
  return fetchJSON(`${BASE}/genes/search?q=${encodeURIComponent(q)}&limit=${limit}`);
}

export async function searchFeatures(q: string, limit = 20): Promise<{ feature: string; count: number; groups: string[] }> {
  return fetchJSON(`${BASE}/features/search?q=${encodeURIComponent(q)}&limit=${limit}`);
}

export async function fetchSharedFeatures(
  groupA: string,
  groupB: string
): Promise<{ group_a: string; group_b: string; count: number; features: string[] }> {
  return fetchJSON(
    `${BASE}/shared-features?group_a=${encodeURIComponent(groupA)}&group_b=${encodeURIComponent(groupB)}`
  );
}

export async function fetchGroupGenes(
  groupName: string,
  limit = 50
): Promise<{ group: string; count: number; genes: string[] }> {
  return fetchJSON(
    `${BASE}/groups/${encodeURIComponent(groupName)}/genes?limit=${limit}`
  );
}

export async function executeSql(query: string): Promise<QueryResult> {
  return fetchJSON(`${BASE}/sql`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ query }),
  });
}

export async function executeCypher(
  query: string,
  metric?: string,
  cutoff?: number
): Promise<QueryResult> {
  return fetchJSON(`${BASE}/cypher`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ query, metric, cutoff }),
  });
}

export async function fetchPairs(
  metric: string,
  cutoff: number,
  limit = 1000,
  offset = 0
): Promise<{ count: number; pairs: Record<string, unknown>[] }> {
  return fetchJSON(
    `${BASE}/pairs?metric=${metric}&cutoff=${cutoff}&limit=${limit}&offset=${offset}`
  );
}

export async function fetchTopPairs(
  metric: string,
  n = 20
): Promise<{ pairs: Record<string, unknown>[] }> {
  return fetchJSON(`${BASE}/top?metric=${metric}&n=${n}`);
}

export async function fetchStatistics(): Promise<Record<string, unknown>> {
  return fetchJSON(`${BASE}/statistics`);
}

export async function fetchMetricSummary(
  metric: string
): Promise<Record<string, number>> {
  return fetchJSON(`${BASE}/statistics/${metric}`);
}

export interface HealthStatus {
  status: "healthy" | "unhealthy";
  timestamp: number;
  version: string;
  details: Record<string, unknown>;
}

export async function fetchHealth(): Promise<HealthStatus> {
  return fetchJSON(`${BASE}/health`);
}

// ── Clustering API ────────────────────────────────────────────────

export interface AlgorithmParameter {
  name: string;
  type: "float" | "int" | "str" | "bool";
  default: number | string | boolean;
  description: string;
  min_value?: number;
  max_value?: number;
  choices?: string[];
}

export interface AlgorithmInfo {
  name: string;
  display_name: string;
  description: string;
  parameters: AlgorithmParameter[];
}

export interface ClusteringResult {
  membership: Record<string, number>;
  num_clusters: number;
  cluster_sizes: Record<string, number>;
  modularity: number | null;
  algorithm: string;
  parameters: Record<string, unknown>;
  metric: string;
  cutoff: number;
}

export async function fetchClusteringAlgorithms(): Promise<{ algorithms: AlgorithmInfo[] }> {
  return fetchJSON(`${BASE}/algorithms/clustering`);
}

export async function runClustering(
  algorithm: string,
  parameters: Record<string, unknown> = {},
  metric?: string,
  cutoff?: number
): Promise<ClusteringResult> {
  return fetchJSON(`${BASE}/graph/cluster`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ algorithm, parameters, metric, cutoff }),
  });
}

export interface ComponentsResult {
  num_components: number;
  components: Record<string, string[]>;
  membership: Record<string, number>;
  min_size_filter: number;
}

export async function fetchComponents(
  metric?: string,
  cutoff?: number,
  minSize: number = 1
): Promise<ComponentsResult> {
  const params = new URLSearchParams();
  if (metric) params.set("metric", metric);
  if (cutoff !== undefined) params.set("cutoff", String(cutoff));
  if (minSize > 1) params.set("min_size", String(minSize));
  return fetchJSON(`${BASE}/graph/components?${params}`);
}

// ── Advanced Filtering API ────────────────────────────────────────

export interface FilterConditionRequest {
  metric: string;
  operator: string;
  value: number | [number, number];
}

export interface FilterResult {
  count: number;
  filters: FilterConditionRequest[];
  logic: string;
  limit: number;
  offset: number;
  pairs: Record<string, unknown>[];
}

export async function filterPairs(
  filters: FilterConditionRequest[],
  logic: "AND" | "OR" = "AND",
  limit: number = 5000,
  offset: number = 0
): Promise<FilterResult> {
  return fetchJSON(`${BASE}/pairs/filter`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ filters, logic, limit, offset }),
  });
}

// ── Export API ────────────────────────────────────────────────────

export type DataExportFormat = "csv" | "tsv" | "json";
export type GraphExportFormat = "graphml" | "gexf" | "json" | "cytoscape";

export function getDataExportUrl(
  format: DataExportFormat,
  metric: string,
  cutoff: number,
  limit: number = 100000
): string {
  return `${BASE}/export/data/${format}?metric=${metric}&cutoff=${cutoff}&limit=${limit}`;
}

export function getGraphExportUrl(
  format: GraphExportFormat,
  metric: string,
  cutoff: number,
  includeLayout: boolean = false,
  includeCommunities: boolean = true
): string {
  const params = new URLSearchParams({
    metric,
    cutoff: String(cutoff),
    include_layout: String(includeLayout),
    include_communities: String(includeCommunities),
  });
  return `${BASE}/export/graph/${format}?${params}`;
}

// ── Biological Analysis API ────────────────────────────────────────

export interface HubGenesResult {
  group: string;
  method: string;
  hops: number;
  genes: Array<{
    gene: string;
    score: number;
    num_diseases?: number;
  }>;
}

export async function fetchHubGenes(
  groupName: string,
  method: "hypergraph" | "edge_weighted" | "projection" = "hypergraph",
  hops: number = 2,
  topN: number = 30
): Promise<HubGenesResult> {
  return fetchJSON(`${BASE}/genes/hub-genes`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      group_name: groupName,
      method,
      hops,
      top_n: topN,
    }),
  });
}

export interface ExplainPairResult {
  group_a: string;
  group_b: string;
  method: string;
  gene_count: number;
  genes: Array<{
    gene: string;
    score: number;
    global_freq?: number;
    specificity?: number;
  }>;
}

export async function explainPair(
  groupA: string,
  groupB: string,
  method: "hypergraph" | "edge_weighted" | "projection" = "hypergraph"
): Promise<ExplainPairResult> {
  return fetchJSON(`${BASE}/genes/explain-pair`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      group_a: groupA,
      group_b: groupB,
      method,
    }),
  });
}

export interface GeneStatistics {
  gene: string;
  group_count: number;
  total_groups: number;
  prevalence_percent: number;
}

export async function fetchGeneStatistics(geneName: string): Promise<GeneStatistics> {
  return fetchJSON(`${BASE}/genes/${encodeURIComponent(geneName)}/statistics`);
}

export interface GeneGroupsResult {
  gene: string;
  group_count: number;
  groups: Array<{ id: number; name: string }>;
}

export async function fetchGeneGroups(
  geneName: string,
  limit: number = 50
): Promise<GeneGroupsResult> {
  return fetchJSON(`${BASE}/genes/${encodeURIComponent(geneName)}/groups?limit=${limit}`);
}

export interface PathResult {
  source: string;
  target: string;
  path_length: number;
  path_nodes: string[];
  connected: boolean;
  shared_genes_along_path?: Array<{
    from: string;
    to: string;
    shared_count: number;
    genes: string[];
  }>;
}

export async function fetchShortestPath(
  source: string,
  target: string,
  metric?: string,
  cutoff?: number
): Promise<PathResult> {
  const params = new URLSearchParams({
    source,
    target,
  });
  if (metric) params.set("metric", metric);
  if (cutoff !== undefined) params.set("cutoff", String(cutoff));
  return fetchJSON(`${BASE}/graph/shortest-path?${params}`);
}

export interface ClusterGenesResult {
  cluster_size: number;
  total_genes: number;
  method: string;
  genes: Array<{
    gene: string;
    score: number;
    in_cluster_count: number;
    cluster_size: number;
    global_count: number;
  }>;
}

export async function analyzeClusterGenes(
  nodeNames: string[],
  method: string = "hypergraph",
  topN: number = 50
): Promise<ClusterGenesResult> {
  return fetchJSON(`${BASE}/genes/cluster-analysis`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      node_names: nodeNames,
      method,
      top_n: topN,
    }),
  });
}

// ── Metric Profile API ────────────────────────────────────────────

export interface MetricProfileResult {
  group: string;
  metrics: Array<{
    metric: string;
    avg: number;
    max: number;
    count: number;
  }>;
}

export async function fetchMetricProfile(
  groupName: string
): Promise<MetricProfileResult> {
  return fetchJSON(`${BASE}/groups/${encodeURIComponent(groupName)}/metric-profile`);
}

/**
 * Helper to check if an error is an APIRequestError
 */
export function isAPIError(error: unknown): error is APIRequestError {
  return error instanceof APIRequestError;
}

/**
 * Helper to safely extract error message from any error
 */
export function getErrorMessage(error: unknown): string {
  if (error instanceof APIRequestError) {
    return error.message;
  }
  if (error instanceof Error) {
    return error.message;
  }
  return String(error);
}

/**
 * Helper to safely extract APIError from any error
 */
export function getAPIError(error: unknown): APIError | null {
  if (error instanceof APIRequestError) {
    return error.apiError;
  }
  return null;
}
