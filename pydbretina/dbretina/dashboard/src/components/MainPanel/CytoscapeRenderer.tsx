import { useRef, useEffect, useCallback, useState } from "react";
import cytoscape, { Core } from "cytoscape";
import fcose from "cytoscape-fcose";
import popper from "cytoscape-popper";
import tippy from "tippy.js";
import "tippy.js/dist/tippy.css";
import { useDashboard } from "../../state/context";
import { GraphNode, LayoutAlgorithm } from "../../state/types";
import { fetchLayout } from "../../api";

cytoscape.use(fcose);
cytoscape.use(popper);

// Map our layout algorithm names to Cytoscape layout names
const LAYOUT_MAP: Record<LayoutAlgorithm, string> = {
  fr: "fcose",
  drl: "fcose", // DRL not available in Cytoscape, use fcose
  kk: "fcose", // KK uses backend positions
  circle: "circle",
  grid: "grid",
};

const COMMUNITY_COLORS = [
  "#4361ee", "#f72585", "#4cc9f0", "#7209b7", "#f77f00",
  "#06d6a0", "#ef476f", "#118ab2", "#ffd166", "#073b4c",
  "#8338ec", "#ff6b6b", "#48bfe3", "#9b5de5", "#00f5d4",
];

export default function CytoscapeRenderer() {
  const containerRef = useRef<HTMLDivElement>(null);
  const cyRef = useRef<Core | null>(null);
  const { state, dispatch } = useDashboard();
  const graphData = state.graphData;
  const [isLayouting, setIsLayouting] = useState(false);

  // Use ref to access current selection state in event handlers (avoid stale closure)
  const selectionRef = useRef(state.selection);
  selectionRef.current = state.selection;

  const handleNodeClick = useCallback(
    (node: GraphNode) => {
      dispatch({ type: "SELECT_NODE", node });
    },
    [dispatch]
  );

  const handleEdgeClick = useCallback(
    (source: string, target: string) => {
      dispatch({ type: "SELECT_EDGE", edge: { source, target } });
    },
    [dispatch]
  );

  // Ref to track if initial layout has been applied
  const initialLayoutApplied = useRef(false);

  // Fetch layout positions from backend for algorithms that need it
  const applyLayout = useCallback(async (
    cy: Core,
    layoutAlgorithm: LayoutAlgorithm,
    metric: string,
    cutoff: number
  ) => {
    const layoutName = LAYOUT_MAP[layoutAlgorithm];

    // For simple layouts, use Cytoscape's built-in
    if (layoutName === "circle" || layoutName === "grid") {
      setIsLayouting(true);
      const layout = cy.layout({ name: layoutName, animate: true, animationDuration: 500 });
      layout.run();
      layout.on("layoutstop", () => setIsLayouting(false));
      return;
    }

    // For force-directed layouts, fetch from backend for better quality
    try {
      setIsLayouting(true);
      const result = await fetchLayout(metric, cutoff, layoutAlgorithm);

      // Apply positions using preset layout
      const positions: Record<string, { x: number; y: number }> = {};
      for (const [nodeId, [x, y]] of Object.entries(result.positions)) {
        positions[nodeId] = { x: x * 500, y: y * 500 }; // Scale positions
      }

      const layout = cy.layout({
        name: "preset",
        positions,
        animate: true,
        animationDuration: 500,
      } as cytoscape.LayoutOptions);
      layout.run();
      layout.on("layoutstop", () => setIsLayouting(false));

      // Store positions in state
      dispatch({ type: "SET_LAYOUT_POSITIONS", positions: result.positions });
    } catch (error) {
      console.error("Failed to fetch layout:", error);
      // Fallback to fcose
      const layout = cy.layout({
        name: "fcose",
        animate: true,
        quality: "default",
        randomize: true,
      } as any);
      layout.run();
      layout.on("layoutstop", () => setIsLayouting(false));
    }
  }, [dispatch]);

  useEffect(() => {
    if (!containerRef.current || !graphData) return;

    // Debug: log container dimensions
    const rect = containerRef.current.getBoundingClientRect();
    console.log("[CytoscapeRenderer] Container dimensions:", rect.width, "x", rect.height);
    console.log("[CytoscapeRenderer] Graph data:", graphData.nodes.length, "nodes,", graphData.edges.length, "edges");

    // Ensure container has dimensions
    if (rect.width === 0 || rect.height === 0) {
      console.warn("[CytoscapeRenderer] Container has zero dimensions, delaying initialization");
      return;
    }

    const elements: cytoscape.ElementDefinition[] = [];

    // Compute max degree for sizing
    const maxDegree = Math.max(...graphData.nodes.map((n) => n.degree), 1);

    for (const node of graphData.nodes) {
      elements.push({
        data: {
          id: node.id,
          label: node.label,
          degree: node.degree,
          community: node.community,
          pagerank: node.pagerank,
          color: COMMUNITY_COLORS[node.community % COMMUNITY_COLORS.length],
          size: 8 + (node.degree / maxDegree) * 32,
        },
      });
    }

    for (const edge of graphData.edges) {
      elements.push({
        data: {
          id: `${edge.source}-${edge.target}`,
          source: edge.source,
          target: edge.target,
          weight: edge.weight,
          shared_features: edge.shared_features,
        },
      });
    }

    const cy = cytoscape({
      container: containerRef.current,
      elements,
      style: [
        {
          selector: "node",
          style: {
            "background-color": "data(color)",
            width: "data(size)",
            height: "data(size)",
            label: "data(label)",
            "font-size": "11px",
            "text-valign": "top",
            "text-halign": "center",
            "text-margin-y": -4,
            "min-zoomed-font-size": 10,
            "border-width": 1,
            "border-color": "#fff",
          },
        },
        {
          selector: "edge",
          style: {
            width: 1.5,
            "line-color": "#ccc",
            opacity: 0.6,
            "curve-style": "bezier",
          },
        },
        {
          selector: "node:selected",
          style: {
            "border-width": 3,
            "border-color": "#ff0",
          },
        },
        {
          selector: ".highlighted",
          style: {
            "border-width": 2,
            "border-color": "#ff6b6b",
          },
        },
        {
          selector: ".dimmed",
          style: {
            opacity: 0.15,
          },
        },
        {
          selector: ".multiselected",
          style: {
            "border-width": 3,
            "border-color": "#4361ee",
            "background-opacity": 0.9,
          },
        },
        {
          selector: ".path-node",
          style: {
            "border-width": 4,
            "border-color": "#f72585",
            "background-opacity": 1,
            "z-index": 999,
          },
        },
        {
          selector: ".path-edge",
          style: {
            "line-color": "#f72585",
            width: 3,
            opacity: 1,
            "z-index": 998,
          },
        },
      ],
      // Use preset layout with random positions initially, then run layout after
      layout: {
        name: "preset",
        positions: () => ({ x: Math.random() * 500, y: Math.random() * 500 }),
      } as cytoscape.LayoutOptions,
      minZoom: 0.05,
      maxZoom: 5,
    });

    cyRef.current = cy;

    cy.on("tap", "node", (evt) => {
      const data = evt.target.data();
      const nodeId = data.id;
      const currentSelection = selectionRef.current;

      // In click mode, toggle selection
      if (currentSelection.mode === "click") {
        // Check if shift key is held for multi-select
        if (evt.originalEvent.shiftKey) {
          if (currentSelection.selectedNodes.has(nodeId)) {
            dispatch({ type: "REMOVE_SELECTED_NODES", nodes: [nodeId] });
          } else {
            dispatch({ type: "ADD_SELECTED_NODES", nodes: [nodeId] });
          }
        } else {
          // Single click without shift - select only this node
          dispatch({ type: "SET_SELECTED_NODES", nodes: new Set([nodeId]) });
        }
      }

      // Also update detail panel
      handleNodeClick({
        id: data.id,
        label: data.label,
        degree: data.degree,
        community: data.community,
        pagerank: data.pagerank,
      });
    });

    cy.on("tap", "edge", (evt) => {
      const data = evt.target.data();
      handleEdgeClick(data.source, data.target);
    });

    // Dismiss selection on background click
    cy.on("tap", (evt) => {
      if (evt.target === cy) {
        dispatch({ type: "SELECT_NODE", node: null });
        dispatch({ type: "SELECT_EDGE", edge: null });
      }
    });

    // Create tooltips for nodes
    cy.nodes().forEach((node) => {
      const ref = (node as any).popperRef();
      const tip = tippy(document.createElement("div"), {
        getReferenceClientRect: ref.getBoundingClientRect,
        trigger: "manual",
        placement: "top",
        content: `
          <strong>${node.data("label")}</strong><br/>
          Degree: ${node.data("degree")}<br/>
          Community: ${node.data("community")}
        `,
        allowHTML: true,
      });

      node.on("mouseover", () => tip.show());
      node.on("mouseout", () => tip.hide());
    });

    // Apply initial layout
    initialLayoutApplied.current = false;
    applyLayout(cy, state.layoutAlgorithm, state.metric, state.cutoff).then(() => {
      initialLayoutApplied.current = true;
    });

    return () => {
      cy.destroy();
    };
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [graphData, handleNodeClick, handleEdgeClick, dispatch]);

  // Handle layout algorithm changes
  useEffect(() => {
    if (cyRef.current && graphData && initialLayoutApplied.current) {
      applyLayout(cyRef.current, state.layoutAlgorithm, state.metric, state.cutoff);
    }
  }, [state.layoutAlgorithm, state.metric, state.cutoff, graphData, applyLayout]);

  // Handle highlight changes
  useEffect(() => {
    const cy = cyRef.current;
    if (!cy) return;

    cy.elements().removeClass("highlighted dimmed");

    if (state.highlightNodes.size > 0) {
      cy.nodes().forEach((n) => {
        if (state.highlightNodes.has(n.id())) {
          n.addClass("highlighted");
        } else {
          n.addClass("dimmed");
        }
      });
      cy.edges().forEach((e) => {
        if (
          !state.highlightNodes.has(e.data("source")) ||
          !state.highlightNodes.has(e.data("target"))
        ) {
          e.addClass("dimmed");
        }
      });
    }
  }, [state.highlightNodes]);

  // Handle selection mode changes
  useEffect(() => {
    const cy = cyRef.current;
    if (!cy) return;

    const mode = state.selection.mode;

    // Enable/disable box selection
    cy.boxSelectionEnabled(mode === "box");

    // Set selection type
    if (mode === "box") {
      cy.selectionType("additive");
    } else {
      cy.selectionType("single");
    }
  }, [state.selection.mode]);

  // Sync selected nodes visual state
  useEffect(() => {
    const cy = cyRef.current;
    if (!cy) return;

    // Update visual selection
    cy.nodes().removeClass("multiselected");
    state.selection.selectedNodes.forEach((nodeId) => {
      const node = cy.getElementById(nodeId);
      if (node.length > 0) {
        node.addClass("multiselected");
      }
    });
  }, [state.selection.selectedNodes]);

  // Handle Cytoscape selection events (box selection)
  useEffect(() => {
    const cy = cyRef.current;
    if (!cy) return;

    const handleBoxEnd = () => {
      if (state.selection.mode !== "box") return;

      const selectedNodeIds = cy.nodes(":selected").map((n) => n.id());
      if (selectedNodeIds.length > 0) {
        dispatch({ type: "ADD_SELECTED_NODES", nodes: selectedNodeIds });
        // Clear Cytoscape's selection (we track in state)
        cy.nodes().unselect();
      }
    };

    cy.on("boxend", handleBoxEnd);

    return () => {
      cy.off("boxend", handleBoxEnd);
    };
  }, [state.selection.mode, dispatch]);

  // Handle path highlighting
  useEffect(() => {
    const cy = cyRef.current;
    if (!cy || !graphData) return;

    // Clear previous path highlighting
    cy.elements().removeClass("path-node path-edge");

    if (state.pathResult?.connected && state.pathResult.path_nodes.length > 0) {
      const pathNodes = state.pathResult.path_nodes.map((n) => n.toLowerCase());

      // Create a map from label to node id
      const labelToId: Record<string, string> = {};
      graphData.nodes.forEach((n) => {
        labelToId[n.label.toLowerCase()] = n.id;
      });

      // Highlight path nodes
      const pathNodeIds: string[] = [];
      pathNodes.forEach((label) => {
        const nodeId = labelToId[label];
        if (nodeId) {
          pathNodeIds.push(nodeId);
          const node = cy.getElementById(nodeId);
          if (node.length > 0) {
            node.addClass("path-node");
          }
        }
      });

      // Highlight path edges (between consecutive path nodes)
      for (let i = 0; i < pathNodeIds.length - 1; i++) {
        const srcId = pathNodeIds[i];
        const tgtId = pathNodeIds[i + 1];
        // Try both directions since edges might be stored either way
        const edge = cy.edges(`[source="${srcId}"][target="${tgtId}"], [source="${tgtId}"][target="${srcId}"]`);
        if (edge.length > 0) {
          edge.addClass("path-edge");
        }
      }

      // Fit view to path if more than 2 nodes
      if (pathNodeIds.length > 2) {
        const pathCollection = cy.collection();
        pathNodeIds.forEach((id) => {
          pathCollection.merge(cy.getElementById(id));
        });
        cy.fit(pathCollection, 50);
      }
    }
  }, [state.pathResult, graphData]);

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
