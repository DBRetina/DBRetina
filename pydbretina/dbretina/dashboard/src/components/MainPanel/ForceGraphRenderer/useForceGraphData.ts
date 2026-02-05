import { useMemo } from "react";
import { DashboardState } from "../../../state/types";
import { ForceGraphData, ForceNode, ForceLink } from "./types";
import { COMMUNITY_COLORS } from "./constants";

/**
 * Hook to transform dashboard GraphData into react-force-graph format.
 * Handles path highlighting and visual properties.
 *
 * NOTE: Selection and highlight state are NOT included here to avoid
 * recomputing graph data on every click. Instead, nodeCanvasObject in
 * the renderer reads selection/highlight from refs for per-frame checks.
 */
export function useForceGraphData(state: DashboardState): ForceGraphData | null {
  return useMemo(() => {
    const { graphData, pathResult, nodeFilter } = state;
    if (!graphData) return null;

    // Apply node filter — when set, only show these nodes and edges between them
    const filteredNodes = nodeFilter
      ? graphData.nodes.filter((n) => nodeFilter.has(n.id))
      : graphData.nodes;
    const filteredEdges = nodeFilter
      ? graphData.edges.filter((e) => nodeFilter.has(e.source) && nodeFilter.has(e.target))
      : graphData.edges;

    // Build path lookup sets for efficient checking
    const pathNodeIds = new Set<string>();
    const pathEdgeKeys = new Set<string>();

    if (pathResult?.connected && pathResult.path_nodes.length > 0) {
      // Build label-to-id map (case-insensitive)
      const labelToId: Record<string, string> = {};
      filteredNodes.forEach((n) => {
        labelToId[n.label.toLowerCase()] = n.id;
      });

      // Collect path node IDs
      pathResult.path_nodes.forEach((label) => {
        const id = labelToId[label.toLowerCase()];
        if (id) pathNodeIds.add(id);
      });

      // Collect path edge keys (both directions)
      const pathIds = Array.from(pathNodeIds);
      for (let i = 0; i < pathIds.length - 1; i++) {
        pathEdgeKeys.add(`${pathIds[i]}-${pathIds[i + 1]}`);
        pathEdgeKeys.add(`${pathIds[i + 1]}-${pathIds[i]}`);
      }
    }

    // Compute max degree for node size normalization
    const maxDegree = Math.max(...filteredNodes.map((n) => n.degree), 1);

    // Compute max shared_features for edge width normalization
    const maxSharedFeatures = Math.max(...filteredEdges.map((e) => e.shared_features), 1);

    // Transform nodes
    // isSelected and isHighlighted are always false here — the renderer
    // checks the latest selection/highlight state from refs each frame.
    const nodes: ForceNode[] = filteredNodes.map((node) => ({
      id: node.id,
      label: node.label,
      degree: node.degree,
      community: node.community,
      pagerank: node.pagerank,
      // Visual properties
      color: COMMUNITY_COLORS[node.community % COMMUNITY_COLORS.length],
      size: 4 + (node.degree / maxDegree) * 20, // Size ranges from 4 to 24
      // State flags — selection/highlight checked at render time via refs
      isHighlighted: false,
      isSelected: false,
      isOnPath: pathNodeIds.has(node.id),
    }));

    // Transform edges to links (react-force-graph uses "links" not "edges")
    const links: ForceLink[] = filteredEdges.map((edge) => ({
      source: edge.source,
      target: edge.target,
      weight: edge.weight,
      shared_features: edge.shared_features,
      // Visual properties - width scales with shared_features (1-5 range)
      color: "rgba(150,150,150,0.4)",
      width: 1 + (edge.shared_features / maxSharedFeatures) * 4,
      isOnPath: pathEdgeKeys.has(`${edge.source}-${edge.target}`),
    }));

    return { nodes, links };
  }, [state.graphData, state.pathResult, state.nodeFilter]);
}
