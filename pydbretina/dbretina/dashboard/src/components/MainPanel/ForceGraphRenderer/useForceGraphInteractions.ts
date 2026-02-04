import { useCallback, MutableRefObject, Dispatch } from "react";
import { DashboardAction, SelectionState, GraphData } from "../../../state/types";
import { ForceNode, ForceLink } from "./types";

interface UseForceGraphInteractionsProps {
  graphRef: MutableRefObject<any>;
  dispatch: Dispatch<DashboardAction>;
  graphData: GraphData | null;
  selectionRef: MutableRefObject<SelectionState>;
}

interface UseForceGraphInteractionsReturn {
  handleNodeClick: (node: ForceNode, event: MouseEvent) => void;
  handleLinkClick: (link: ForceLink, event: MouseEvent) => void;
  handleBackgroundClick: () => void;
}

/**
 * Hook to handle user interactions with the force graph.
 * Provides handlers for node clicks, link clicks, and background clicks.
 */
export function useForceGraphInteractions({
  dispatch,
  selectionRef,
}: UseForceGraphInteractionsProps): UseForceGraphInteractionsReturn {
  /**
   * Handle node click - supports single select and multi-select (shift+click)
   */
  const handleNodeClick = useCallback(
    (node: ForceNode, event: MouseEvent) => {
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
    [dispatch, selectionRef]
  );

  /**
   * Handle link/edge click - shows edge detail panel
   */
  const handleLinkClick = useCallback(
    (link: ForceLink, _event: MouseEvent) => {
      // Extract source/target IDs (links may have resolved node objects)
      const source = typeof link.source === "string" ? link.source : link.source.id;
      const target = typeof link.target === "string" ? link.target : link.target.id;

      dispatch({ type: "SELECT_EDGE", edge: { source, target } });
    },
    [dispatch]
  );

  /**
   * Handle background click - clears selection
   */
  const handleBackgroundClick = useCallback(() => {
    dispatch({ type: "SELECT_NODE", node: null });
    dispatch({ type: "SELECT_EDGE", edge: null });
  }, [dispatch]);

  return {
    handleNodeClick,
    handleLinkClick,
    handleBackgroundClick,
  };
}
