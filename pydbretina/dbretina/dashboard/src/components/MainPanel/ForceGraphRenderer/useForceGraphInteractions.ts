import { useCallback, useRef, MutableRefObject, Dispatch } from "react";
import { DashboardAction, SelectionState, CompareState, GraphData } from "../../../state/types";
import { ForceNode, ForceLink } from "./types";

interface UseForceGraphInteractionsProps {
  graphRef: MutableRefObject<any>;
  dispatch: Dispatch<DashboardAction>;
  graphData: GraphData | null;
  selectionRef: MutableRefObject<SelectionState>;
  compareStateRef: MutableRefObject<CompareState>;
}

interface UseForceGraphInteractionsReturn {
  handleNodeClick: (node: ForceNode, event: MouseEvent) => void;
  handleLinkClick: (link: ForceLink, event: MouseEvent) => void;
  handleBackgroundClick: () => void;
}

/**
 * Hook to handle user interactions with the force graph.
 * Provides handlers for node clicks, link clicks, and background clicks.
 *
 * Uses a timestamp guard to prevent background clicks from immediately
 * clearing a selection set by a node/link click (race condition in
 * react-force-graph with dense graphs).
 */
export function useForceGraphInteractions({
  dispatch,
  selectionRef,
  compareStateRef,
}: UseForceGraphInteractionsProps): UseForceGraphInteractionsReturn {
  // Guard: track last element click time to avoid background click race condition
  const lastElementClickTime = useRef(0);

  /**
   * Handle node click - supports single select, multi-select (shift+click),
   * and compare mode (picks second node).
   */
  const handleNodeClick = useCallback(
    (node: ForceNode, event: MouseEvent) => {
      lastElementClickTime.current = Date.now();

      // Compare mode: if waiting for nodeB, set it instead of normal selection
      const cmp = compareStateRef.current;
      if (cmp.active && !cmp.nodeB) {
        dispatch({
          type: "SET_COMPARE_NODE_B",
          nodeB: {
            id: node.id,
            label: node.label,
            degree: node.degree,
            community: node.community,
            pagerank: node.pagerank,
          },
        });
        return;
      }

      const currentSelection = selectionRef.current;

      if (event.shiftKey) {
        // Multi-select: toggle this node
        if (currentSelection.selectedNodes.has(node.id)) {
          dispatch({ type: "REMOVE_SELECTED_NODES", nodes: [node.id] });
        } else {
          dispatch({ type: "ADD_SELECTED_NODES", nodes: [node.id] });
        }
      } else {
        // Single select: select only this node
        dispatch({ type: "SET_SELECTED_NODES", nodes: new Set([node.id]) });
      }

      // Update detail panel with node info
      dispatch({
        type: "SELECT_NODE",
        node: {
          id: node.id,
          label: node.label,
          degree: node.degree,
          community: node.community,
          pagerank: node.pagerank,
        },
      });
    },
    [dispatch, selectionRef, compareStateRef]
  );

  /**
   * Handle link/edge click - shows edge detail panel
   */
  const handleLinkClick = useCallback(
    (link: ForceLink, _event: MouseEvent) => {
      lastElementClickTime.current = Date.now();
      // Extract source/target IDs (links may have resolved node objects)
      const source = typeof link.source === "string" ? link.source : link.source.id;
      const target = typeof link.target === "string" ? link.target : link.target.id;

      dispatch({ type: "SELECT_EDGE", edge: { source, target } });
    },
    [dispatch]
  );

  /**
   * Handle background click - clears selection.
   * Ignores clicks within 200ms of a node/link click to prevent race conditions.
   */
  const handleBackgroundClick = useCallback(() => {
    if (Date.now() - lastElementClickTime.current < 200) return;
    dispatch({ type: "SELECT_NODE", node: null });
    dispatch({ type: "SELECT_EDGE", edge: null });
  }, [dispatch]);

  return {
    handleNodeClick,
    handleLinkClick,
    handleBackgroundClick,
  };
}
