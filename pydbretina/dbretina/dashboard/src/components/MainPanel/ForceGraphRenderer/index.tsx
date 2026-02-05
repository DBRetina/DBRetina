import { useRef, useState, useEffect, useCallback } from "react";
import ForceGraph2D from "react-force-graph-2d";
import ForceGraph3D from "react-force-graph-3d";
import { useDashboard } from "../../../state/context";
import { useForceGraphData } from "./useForceGraphData";
import { useForceGraphInteractions } from "./useForceGraphInteractions";
import { ForceNode, ForceLink } from "./types";
import ForceGraphControls from "./ForceGraphControls";
import { THRESHOLDS, COLORS, PHYSICS } from "./constants";
import { fetchLayout } from "../../../api";
import Graph from "graphology";
import forceAtlas2 from "graphology-layout-forceatlas2";

/**
 * ForceGraphRenderer - High-performance graph visualization using react-force-graph.
 *
 * Automatically switches between 2D (Canvas) and 3D (WebGL) based on graph size:
 * - < 3000 nodes: 2D Canvas (better for interaction, labels)
 * - >= 3000 nodes: 3D WebGL (better for performance)
 */
export default function ForceGraphRenderer() {
  const containerRef = useRef<HTMLDivElement>(null);
  const graphRef = useRef<any>();
  const { state, dispatch } = useDashboard();

  // Track selection state in ref to avoid stale closures in event handlers
  const selectionRef = useRef(state.selection);
  selectionRef.current = state.selection;

  // Track compare state in ref for interaction handlers
  const compareStateRef = useRef(state.compareState);
  compareStateRef.current = state.compareState;

  // Track highlight nodes in ref for per-frame rendering
  const highlightNodesRef = useRef(state.highlightNodes);
  highlightNodesRef.current = state.highlightNodes;

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
    compareStateRef,
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

  // Configure d3-force physics when graph instance is ready
  useEffect(() => {
    const graph = graphRef.current;
    if (!graph || !forceGraphData) return;

    const nodeCount = forceGraphData.nodes.length;
    const isLarge = nodeCount > THRESHOLDS.REDUCE_PHYSICS;

    // Configure charge (repulsion) force
    const charge = graph.d3Force("charge");
    if (charge) {
      charge
        .strength(isLarge ? -30 : PHYSICS.CHARGE_STRENGTH)
        .distanceMax(isLarge ? 300 : PHYSICS.CHARGE_DISTANCE_MAX);
    }

    // Configure link (spring) force
    const link = graph.d3Force("link");
    if (link) {
      link
        .distance(isLarge ? 30 : PHYSICS.LINK_DISTANCE)
        .strength(PHYSICS.LINK_STRENGTH);
    }

    // Configure center force
    const center = graph.d3Force("center");
    if (center) {
      center.strength(PHYSICS.CENTER_STRENGTH);
    }

    graph.d3ReheatSimulation?.();
  }, [forceGraphData]);

  // Auto-fit graph when data changes
  useEffect(() => {
    if (forceGraphData && graphRef.current) {
      const timer = setTimeout(() => {
        graphRef.current?.zoomToFit?.(400, 50);
      }, 800);
      return () => clearTimeout(timer);
    }
  }, [forceGraphData]);

  // Handle layout algorithm changes
  useEffect(() => {
    if (!forceGraphData || !graphRef.current) return;

    const nodes = forceGraphData.nodes;
    const nc = nodes.length;
    const layout = state.layoutAlgorithm;

    // Helper: clear all fixed positions and let simulation run
    const clearFixed = () => {
      nodes.forEach((node) => {
        node.fx = undefined;
        node.fy = undefined;
      });
    };

    // Helper: apply positions and fit
    const applyPositionsAndFit = (
      positions: Record<string, [number, number]>,
      scale: number
    ) => {
      nodes.forEach((node) => {
        const pos = positions[node.id];
        if (pos) {
          node.x = pos[0] * scale;
          node.y = pos[1] * scale;
          node.fx = pos[0] * scale;
          node.fy = pos[1] * scale;
        }
      });
      graphRef.current?.d3ReheatSimulation?.();
      setTimeout(() => {
        // Release fixed positions for interaction
        clearFixed();
        graphRef.current?.zoomToFit?.(400, 50);
      }, 300);
    };

    // "force" = live d3-force simulation (default)
    if (layout === "force") {
      clearFixed();
      graphRef.current?.d3ReheatSimulation?.();
      setTimeout(() => graphRef.current?.zoomToFit?.(400, 50), 800);
      return;
    }

    // Circle layout: arrange nodes in a circle
    if (layout === "circle") {
      const radius = Math.max(100, nc * 3);
      nodes.forEach((node, i) => {
        const angle = (2 * Math.PI * i) / nc;
        node.fx = radius * Math.cos(angle);
        node.fy = radius * Math.sin(angle);
      });
      graphRef.current?.d3ReheatSimulation?.();
      setTimeout(() => graphRef.current?.zoomToFit?.(400, 50), 100);
      return;
    }

    // Grid layout: arrange nodes in a grid
    if (layout === "grid") {
      const cols = Math.ceil(Math.sqrt(nc));
      const spacing = 50;
      nodes.forEach((node, i) => {
        node.fx = (i % cols) * spacing - (cols * spacing) / 2;
        node.fy = Math.floor(i / cols) * spacing - (Math.ceil(nc / cols) * spacing) / 2;
      });
      graphRef.current?.d3ReheatSimulation?.();
      setTimeout(() => graphRef.current?.zoomToFit?.(400, 50), 100);
      return;
    }

    // ForceAtlas2: compute layout client-side using graphology
    if (layout === "fa2") {
      clearFixed();

      // Build graphology graph
      const g = new Graph();
      nodes.forEach((n) => g.addNode(n.id, { x: n.x ?? Math.random() * 500, y: n.y ?? Math.random() * 500 }));
      forceGraphData.links.forEach((l) => {
        const src = typeof l.source === "string" ? l.source : l.source.id;
        const tgt = typeof l.target === "string" ? l.target : l.target.id;
        if (g.hasNode(src) && g.hasNode(tgt) && !g.hasEdge(src, tgt)) {
          g.addEdge(src, tgt, { weight: l.weight });
        }
      });

      // Run ForceAtlas2
      const iterations = nc < 500 ? 200 : nc < 2000 ? 100 : 50;
      forceAtlas2.assign(g, {
        iterations,
        settings: {
          gravity: 1,
          scalingRatio: nc < 200 ? 10 : 5,
          strongGravityMode: true,
          barnesHutOptimize: nc > 500,
          barnesHutTheta: 0.5,
        },
      });

      // Apply positions
      const positions: Record<string, [number, number]> = {};
      g.forEachNode((id, attrs) => {
        positions[id] = [attrs.x, attrs.y];
      });
      applyPositionsAndFit(positions, 1);
      return;
    }

    // Backend-computed layouts: FR, DRL, KK
    if (layout === "fr" || layout === "drl" || layout === "kk") {
      clearFixed();

      fetchLayout(state.metric, state.cutoff, layout)
        .then((result) => {
          if (!graphRef.current) return;
          applyPositionsAndFit(result.positions, 500);
        })
        .catch(() => {
          // Fallback: reheat simulation
          graphRef.current?.d3ReheatSimulation?.();
          setTimeout(() => graphRef.current?.zoomToFit?.(400, 50), 800);
        });
    }
  }, [state.layoutAlgorithm]);

  const nodeCount = state.graphData?.nodes.length ?? 0;
  const use3D = nodeCount >= THRESHOLDS.USE_3D;

  /**
   * Custom 2D node rendering with Canvas API.
   * Reads selection/highlight state from refs for per-frame accuracy
   * without triggering graph data recomputation on clicks.
   */
  const nodeCanvasObject = useCallback(
    (node: ForceNode, ctx: CanvasRenderingContext2D, globalScale: number) => {
      const size = node.size;
      const x = node.x!;
      const y = node.y!;

      // Check selection/highlight state from refs (always current, no data recompute)
      const isSelected = selectionRef.current.selectedNodes.has(node.id);
      const isHighlighted = highlightNodesRef.current.has(node.id);

      // Draw node circle
      ctx.beginPath();
      ctx.arc(x, y, size, 0, 2 * Math.PI);
      ctx.fillStyle = node.isOnPath ? COLORS.PATH_HIGHLIGHT : node.color;
      ctx.fill();

      // Draw border for selection/path highlighting
      if (isSelected) {
        ctx.strokeStyle = COLORS.SELECTION_HIGHLIGHT;
        ctx.lineWidth = 3 / globalScale;
        ctx.stroke();
      } else if (node.isOnPath) {
        ctx.strokeStyle = COLORS.PATH_HIGHLIGHT;
        ctx.lineWidth = 4 / globalScale;
        ctx.stroke();
      } else if (isHighlighted) {
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
      ctx.arc(node.x!, node.y!, node.size + 8, 0, 2 * Math.PI);
      ctx.fillStyle = color;
      ctx.fill();
    },
    []
  );

  // Show empty container while loading
  if (!forceGraphData) {
    return <div className="graph-container" />;
  }

  const isLarge = nodeCount > THRESHOLDS.REDUCE_PHYSICS;

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
    linkHoverPrecision: 6, // Wider hit area for edge clicks
    onNodeClick: interactions.handleNodeClick,
    onNodeHover: setHoveredNode,
    onLinkClick: interactions.handleLinkClick,
    onBackgroundClick: interactions.handleBackgroundClick,
    // Physics settings
    d3AlphaDecay: isLarge ? 0.1 : PHYSICS.ALPHA_DECAY,
    d3VelocityDecay: PHYSICS.VELOCITY_DECAY,
    warmupTicks: isLarge ? 0 : PHYSICS.WARMUP_TICKS,
    cooldownTicks: isLarge ? 0 : PHYSICS.COOLDOWN_TICKS,
  };

  return (
    <div
      ref={containerRef}
      className="graph-container"
      style={{ position: "relative", height: "100%" }}
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
