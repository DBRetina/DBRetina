import { useRef, useEffect, useCallback, useState } from "react";
import Graph from "graphology";
import Sigma from "sigma";
import { useDashboard } from "../../state/context";
import { GraphNode } from "../../state/types";
import { fetchLayout } from "../../api";

const COMMUNITY_COLORS = [
  "#4361ee", "#f72585", "#4cc9f0", "#7209b7", "#f77f00",
  "#06d6a0", "#ef476f", "#118ab2", "#ffd166", "#073b4c",
  "#8338ec", "#ff6b6b", "#48bfe3", "#9b5de5", "#00f5d4",
];

export default function SigmaRenderer() {
  const containerRef = useRef<HTMLDivElement>(null);
  const sigmaRef = useRef<Sigma | null>(null);
  const graphRef = useRef<Graph | null>(null);
  const { state, dispatch } = useDashboard();
  const graphData = state.graphData;
  const [isLayouting, setIsLayouting] = useState(false);

  // Use ref to access current selection state in event handlers
  const selectionRef = useRef(state.selection);
  selectionRef.current = state.selection;

  const handleNodeClick = useCallback(
    (node: GraphNode, shiftKey: boolean = false) => {
      const currentSelection = selectionRef.current;

      // In click mode, handle selection
      if (currentSelection.mode === "click") {
        if (shiftKey) {
          if (currentSelection.selectedNodes.has(node.id)) {
            dispatch({ type: "REMOVE_SELECTED_NODES", nodes: [node.id] });
          } else {
            dispatch({ type: "ADD_SELECTED_NODES", nodes: [node.id] });
          }
        } else {
          dispatch({ type: "SET_SELECTED_NODES", nodes: new Set([node.id]) });
        }
      }

      // Also update detail panel
      dispatch({ type: "SELECT_NODE", node });
    },
    [dispatch]
  );

  // Ref to track if initial layout has been applied
  const initialLayoutApplied = useRef(false);

  // Fetch and apply layout positions
  const applyLayout = useCallback(async (
    graph: Graph,
    layoutAlgorithm: string,
    metric: string,
    cutoff: number
  ) => {
    try {
      setIsLayouting(true);
      const result = await fetchLayout(metric, cutoff, layoutAlgorithm);

      // Apply positions to graph nodes
      for (const [nodeId, [x, y]] of Object.entries(result.positions)) {
        if (graph.hasNode(nodeId)) {
          graph.setNodeAttribute(nodeId, "x", x * 500);
          graph.setNodeAttribute(nodeId, "y", y * 500);
        }
      }

      // Refresh sigma renderer
      sigmaRef.current?.refresh();

      // Store positions in state
      dispatch({ type: "SET_LAYOUT_POSITIONS", positions: result.positions });
    } catch (error) {
      console.error("Failed to fetch layout:", error);
      // Keep random positions as fallback
    } finally {
      setIsLayouting(false);
    }
  }, [dispatch]);

  useEffect(() => {
    if (!containerRef.current || !graphData) return;

    const graph = new Graph();
    const maxDegree = Math.max(...graphData.nodes.map((n) => n.degree), 1);

    for (const node of graphData.nodes) {
      graph.addNode(node.id, {
        label: node.label,
        x: Math.random() * 1000,
        y: Math.random() * 1000,
        size: 2 + (node.degree / maxDegree) * 12,
        color: COMMUNITY_COLORS[node.community % COMMUNITY_COLORS.length],
        degree: node.degree,
        community: node.community,
        pagerank: node.pagerank,
      });
    }

    for (const edge of graphData.edges) {
      const key = `${edge.source}-${edge.target}`;
      if (!graph.hasEdge(key) && graph.hasNode(edge.source) && graph.hasNode(edge.target)) {
        graph.addEdgeWithKey(key, edge.source, edge.target, {
          weight: edge.weight,
          color: "rgba(200,200,200,0.3)",
        });
      }
    }

    const renderer = new Sigma(graph, containerRef.current, {
      renderEdgeLabels: false,
      defaultEdgeType: "line",
      labelRenderedSizeThreshold: 8,
      labelFont: "sans-serif",
      labelSize: 10,
    });

    renderer.on("clickNode", ({ node, event }) => {
      const attrs = graph.getNodeAttributes(node);
      handleNodeClick(
        {
          id: node,
          label: attrs.label,
          degree: attrs.degree,
          community: attrs.community,
          pagerank: attrs.pagerank,
        },
        event.original.shiftKey
      );
    });

    graphRef.current = graph;
    sigmaRef.current = renderer;

    // Fetch and apply initial layout
    initialLayoutApplied.current = false;
    applyLayout(graph, state.layoutAlgorithm, state.metric, state.cutoff).then(() => {
      initialLayoutApplied.current = true;
    });

    return () => {
      renderer.kill();
    };
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [graphData, handleNodeClick]);

  // Handle layout algorithm changes
  useEffect(() => {
    if (graphRef.current && graphData && initialLayoutApplied.current) {
      applyLayout(graphRef.current, state.layoutAlgorithm, state.metric, state.cutoff);
    }
  }, [state.layoutAlgorithm, state.metric, state.cutoff, graphData, applyLayout]);

  // Handle selection and path visual updates
  useEffect(() => {
    const renderer = sigmaRef.current;
    const graph = graphRef.current;
    if (!renderer || !graph) return;

    // Build path node set for quick lookup
    const pathNodeIds = new Set<string>();
    const pathEdges = new Set<string>();

    if (state.pathResult?.connected && state.pathResult.path_nodes.length > 0 && graphData) {
      const pathNodeLabels = state.pathResult.path_nodes.map((n) => n.toLowerCase());

      // Build label to id map
      const labelToId: Record<string, string> = {};
      graphData.nodes.forEach((n) => {
        labelToId[n.label.toLowerCase()] = n.id;
      });

      // Get path node IDs
      pathNodeLabels.forEach((label) => {
        const nodeId = labelToId[label];
        if (nodeId) pathNodeIds.add(nodeId);
      });

      // Get path edge IDs
      const pathNodeIdList = Array.from(pathNodeIds);
      for (let i = 0; i < pathNodeIdList.length - 1; i++) {
        const srcId = pathNodeIdList[i];
        const tgtId = pathNodeIdList[i + 1];
        pathEdges.add(`${srcId}-${tgtId}`);
        pathEdges.add(`${tgtId}-${srcId}`);
      }
    }

    // Update node reducer to show selection and path
    renderer.setSetting("nodeReducer", (node, data) => {
      const isSelected = state.selection.selectedNodes.has(node);
      const isOnPath = pathNodeIds.has(node);

      if (isOnPath) {
        return {
          ...data,
          color: "#f72585", // Path highlight color
          size: data.size * 1.5,
          zIndex: 2,
        };
      }
      if (isSelected) {
        return {
          ...data,
          color: "#4361ee", // Selection highlight color
          zIndex: 1,
        };
      }
      // Dim other nodes if path is active
      if (pathNodeIds.size > 0) {
        return {
          ...data,
          color: data.color,
          // Make nodes slightly dimmer when path is shown
        };
      }
      return data;
    });

    // Update edge reducer to show path edges
    renderer.setSetting("edgeReducer", (edge, data) => {
      if (pathEdges.has(edge)) {
        return {
          ...data,
          color: "#f72585",
          size: 3,
          zIndex: 1,
        };
      }
      // Dim other edges if path is active
      if (pathEdges.size > 0) {
        return {
          ...data,
          color: "rgba(200,200,200,0.1)",
        };
      }
      return data;
    });

    renderer.refresh();
  }, [state.selection.selectedNodes, state.pathResult, graphData]);

  return (
    <div className="graph-container" style={{ position: "relative", width: "100%", height: "100%" }}>
      <div ref={containerRef} style={{ position: "absolute", top: 0, left: 0, right: 0, bottom: 0 }} />
      {isLayouting && (
        <div
          style={{
            position: "absolute",
            top: 8,
            right: 8,
            background: "rgba(0,0,0,0.7)",
            color: "white",
            padding: "4px 12px",
            borderRadius: 4,
            fontSize: 12,
          }}
        >
          Computing layout...
        </div>
      )}
    </div>
  );
}
