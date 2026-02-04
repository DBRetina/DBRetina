import { useRef, useState, useEffect, useCallback } from "react";
import ForceGraph2D from "react-force-graph-2d";
import ForceGraph3D from "react-force-graph-3d";
import { useDashboard } from "../../../state/context";
import { useForceGraphData } from "./useForceGraphData";
import { useForceGraphInteractions } from "./useForceGraphInteractions";
import { ForceNode, ForceLink } from "./types";
import ForceGraphControls from "./ForceGraphControls";
import { THRESHOLDS, COLORS } from "./constants";

/**
 * ForceGraphRenderer - High-performance graph visualization using react-force-graph.
 *
 * Automatically switches between 2D (Canvas) and 3D (WebGL) based on graph size:
 * - < 3000 nodes: 2D Canvas (better for interaction, labels)
 * - >= 3000 nodes: 3D WebGL (better for performance)
 */
export default function ForceGraphRenderer() {
  const containerRef = useRef<HTMLDivElement>(null);
  // Use any for graph ref since the library types are complex with generics
  const graphRef = useRef<any>();
  const { state, dispatch } = useDashboard();

  // Track selection state in ref to avoid stale closures in event handlers
  const selectionRef = useRef(state.selection);
  selectionRef.current = state.selection;

  const [dimensions, setDimensions] = useState({ width: 800, height: 600 });
  const [hoveredNode, setHoveredNode] = useState<ForceNode | null>(null);

  // Transform data for force-graph
  const forceGraphData = useForceGraphData(state);

  // Get interaction handlers
  const interactions = useForceGraphInteractions({
    graphRef,
    dispatch,
    graphData: state.graphData,
    selectionRef,
  });

  // Track container dimensions with ResizeObserver
  useEffect(() => {
    if (!containerRef.current) return;

    const observer = new ResizeObserver((entries) => {
      const { width, height } = entries[0].contentRect;
      if (width > 0 && height > 0) {
        setDimensions({ width, height });
      }
    });

    observer.observe(containerRef.current);
    return () => observer.disconnect();
  }, []);

  // Auto-fit graph when data changes
  useEffect(() => {
    if (forceGraphData && graphRef.current) {
      // Delay to let physics settle a bit first
      const timer = setTimeout(() => {
        graphRef.current?.zoomToFit?.(400, 50);
      }, 500);
      return () => clearTimeout(timer);
    }
  }, [forceGraphData]);

  const nodeCount = state.graphData?.nodes.length ?? 0;
  const use3D = nodeCount >= THRESHOLDS.USE_3D;

  /**
   * Custom 2D node rendering with Canvas API
   */
  const nodeCanvasObject = useCallback(
    (node: ForceNode, ctx: CanvasRenderingContext2D, globalScale: number) => {
      const size = node.size;
      const x = node.x!;
      const y = node.y!;

      // Draw node circle
      ctx.beginPath();
      ctx.arc(x, y, size, 0, 2 * Math.PI);
      ctx.fillStyle = node.isOnPath ? COLORS.PATH_HIGHLIGHT : node.color;
      ctx.fill();

      // Draw border for selection/path highlighting
      if (node.isSelected) {
        ctx.strokeStyle = COLORS.SELECTION_HIGHLIGHT;
        ctx.lineWidth = 3 / globalScale;
        ctx.stroke();
      } else if (node.isOnPath) {
        ctx.strokeStyle = COLORS.PATH_HIGHLIGHT;
        ctx.lineWidth = 4 / globalScale;
        ctx.stroke();
      } else if (node.isHighlighted) {
        ctx.strokeStyle = "#ff6b6b";
        ctx.lineWidth = 2 / globalScale;
        ctx.stroke();
      } else {
        ctx.strokeStyle = COLORS.NODE_BORDER;
        ctx.lineWidth = 1 / globalScale;
        ctx.stroke();
      }

      // Draw label when zoomed in enough (and not too many nodes)
      if (globalScale > 0.5 && nodeCount < THRESHOLDS.DISABLE_LABELS) {
        const fontSize = Math.max(10 / globalScale, 4);
        ctx.font = `${fontSize}px sans-serif`;
        ctx.textAlign = "center";
        ctx.textBaseline = "top";
        ctx.fillStyle = COLORS.TEXT;
        ctx.fillText(node.label, x, y + size + 2);
      }
    },
    [nodeCount]
  );

  /**
   * Custom pointer area for better click detection
   */
  const nodePointerAreaPaint = useCallback(
    (node: ForceNode, color: string, ctx: CanvasRenderingContext2D) => {
      ctx.beginPath();
      ctx.arc(node.x!, node.y!, node.size + 4, 0, 2 * Math.PI);
      ctx.fillStyle = color;
      ctx.fill();
    },
    []
  );

  // Show empty container while loading
  if (!forceGraphData) {
    return <div className="graph-container" />;
  }

  // Common props for both 2D and 3D
  const commonProps = {
    graphData: forceGraphData,
    width: dimensions.width,
    height: dimensions.height,
    nodeId: "id" as const,
    nodeLabel: "label" as const,
    nodeColor: (n: ForceNode) => (n.isOnPath ? COLORS.PATH_HIGHLIGHT : n.color),
    nodeVal: (n: ForceNode) => n.size,
    linkColor: (l: ForceLink) => (l.isOnPath ? COLORS.LINK_PATH : l.color),
    linkWidth: (l: ForceLink) => (l.isOnPath ? 3 : l.width),
    onNodeClick: interactions.handleNodeClick,
    onNodeHover: setHoveredNode,
    onLinkClick: interactions.handleLinkClick,
    onBackgroundClick: interactions.handleBackgroundClick,
    // Physics settings - reduce for large graphs
    d3AlphaDecay: nodeCount > THRESHOLDS.REDUCE_PHYSICS ? 0.1 : 0.02,
    d3VelocityDecay: 0.3,
    warmupTicks: nodeCount > THRESHOLDS.REDUCE_PHYSICS ? 0 : 100,
    cooldownTicks: nodeCount > THRESHOLDS.REDUCE_PHYSICS ? 0 : 200,
  };

  return (
    <div
      ref={containerRef}
      className="graph-container"
      style={{ position: "relative", width: "100%", height: "100%" }}
    >
      {use3D ? (
        <ForceGraph3D
          ref={graphRef as any}
          {...commonProps}
          linkOpacity={0.4}
          nodeOpacity={0.9}
        />
      ) : (
        <ForceGraph2D
          ref={graphRef}
          {...commonProps}
          nodeCanvasObject={nodeCanvasObject}
          nodePointerAreaPaint={nodePointerAreaPaint}
        />
      )}

      <ForceGraphControls graphRef={graphRef} />

      {/* Tooltip for hovered node */}
      {hoveredNode && (
        <div
          style={{
            position: "absolute",
            top: 8,
            left: 8,
            background: "rgba(0,0,0,0.85)",
            color: "white",
            padding: "8px 12px",
            borderRadius: 4,
            fontSize: 12,
            pointerEvents: "none",
            maxWidth: 250,
            zIndex: 100,
          }}
        >
          <strong>{hoveredNode.label}</strong>
          <br />
          Degree: {hoveredNode.degree}
          <br />
          Community: {hoveredNode.community}
          {hoveredNode.pagerank > 0 && (
            <>
              <br />
              PageRank: {hoveredNode.pagerank.toFixed(4)}
            </>
          )}
        </div>
      )}

      {/* Mode indicator */}
      <div
        style={{
          position: "absolute",
          top: 8,
          right: 60,
          background: "rgba(0,0,0,0.5)",
          color: "white",
          padding: "4px 8px",
          borderRadius: 4,
          fontSize: 10,
          pointerEvents: "none",
        }}
      >
        {use3D ? "3D WebGL" : "2D Canvas"}
      </div>
    </div>
  );
}
