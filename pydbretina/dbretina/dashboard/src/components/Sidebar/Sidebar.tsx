import { useCallback } from "react";
import { useDashboard } from "../../state/context";
import { GraphEdge } from "../../state/types";
import DatasetInfo from "./DatasetInfo";
import MetricSelector from "./MetricSelector";
import CutoffSlider from "./CutoffSlider";
import GroupSearch from "./GroupSearch";
import FeatureSearch from "./FeatureSearch";
import GeneExplorer from "./GeneExplorer";
import HubGenesPanel from "./HubGenesPanel";
import PathFinder from "./PathFinder";
import CommunityPanel from "./CommunityPanel";
import ClusteringPanel from "./ClusteringPanel";
import ExportPanel from "./ExportPanel";
import { FilterBuilder, FilterCondition } from "./AdvancedFilter";

/** Test a single edge metric value against one condition. Returns false when
 * the edge lacks the metric (undefined) so missing data never "matches". */
function conditionMatches(edge: GraphEdge, cond: FilterCondition): boolean {
  const v = edge[cond.metric as keyof GraphEdge] as number | undefined;
  if (v === undefined || v === null || Number.isNaN(v)) return false;
  if (cond.operator === "between") {
    const [min, max] = Array.isArray(cond.value) ? cond.value : [cond.value, cond.value];
    return v >= min && v <= max;
  }
  const t = Array.isArray(cond.value) ? cond.value[0] : cond.value;
  switch (cond.operator) {
    case ">=":
      return v >= t;
    case "<=":
      return v <= t;
    case ">":
      return v > t;
    case "<":
      return v < t;
    case "==":
      return v === t;
    case "!=":
      return v !== t;
    default:
      return false;
  }
}

export default function Sidebar() {
  const { state, dispatch } = useDashboard();

  // Apply the Advanced Filter as a client-side edge filter over the loaded
  // graph. An edge matches when (AND) every condition holds / (OR) any holds.
  const handleFilterApply = useCallback(
    (filters: FilterCondition[], logic: "AND" | "OR") => {
      const edges = state.graphData?.edges;
      if (!edges || filters.length === 0) {
        dispatch({ type: "CLEAR_EDGE_FILTER" });
        dispatch({ type: "CLEAR_NODE_FILTER" });
        return;
      }
      const matchEdge = (e: GraphEdge) =>
        logic === "AND"
          ? filters.every((c) => conditionMatches(e, c))
          : filters.some((c) => conditionMatches(e, c));

      const edgeKeys = new Set<string>();
      const nodeKeys = new Set<string>();
      for (const e of edges) {
        if (matchEdge(e)) {
          edgeKeys.add(`${e.source}|${e.target}`);
          nodeKeys.add(e.source);
          nodeKeys.add(e.target);
        }
      }
      dispatch({ type: "SET_EDGE_FILTER", edges: edgeKeys });
      dispatch({ type: "SET_NODE_FILTER", nodes: nodeKeys });
    },
    [state.graphData, dispatch]
  );

  const handleFilterClear = useCallback(() => {
    dispatch({ type: "CLEAR_EDGE_FILTER" });
    dispatch({ type: "CLEAR_NODE_FILTER" });
  }, [dispatch]);

  return (
    <div className="sidebar">
      <div className="sidebar-header">
        <h1>DBRetina</h1>
      </div>
      <DatasetInfo />
      <MetricSelector />
      <CutoffSlider />
      <GroupSearch />
      {state.info?.has_genes && <FeatureSearch />}
      {state.info?.has_genes && (
        <>
          <GeneExplorer />
          <HubGenesPanel />
          <PathFinder />
        </>
      )}
      <FilterBuilder onApply={handleFilterApply} onClear={handleFilterClear} />
      <ClusteringPanel />
      <ExportPanel />
      <CommunityPanel />
    </div>
  );
}
