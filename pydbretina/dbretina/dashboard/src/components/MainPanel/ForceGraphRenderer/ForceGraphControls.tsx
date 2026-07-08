import { MutableRefObject } from "react";

interface Props {
  graphRef: MutableRefObject<any>;
}

/**
 * Control buttons for the force graph (zoom, fit, reheat, etc.)
 */
export default function ForceGraphControls({ graphRef }: Props) {
  const handleZoomIn = () => {
    const graph = graphRef.current;
    if (graph) {
      const currentZoom = graph.zoom?.() ?? 1;
      graph.zoom?.(currentZoom * 1.5, 300);
    }
  };

  const handleZoomOut = () => {
    const graph = graphRef.current;
    if (graph) {
      const currentZoom = graph.zoom?.() ?? 1;
      graph.zoom?.(currentZoom / 1.5, 300);
    }
  };

  const handleFit = () => {
    graphRef.current?.zoomToFit?.(400, 50);
  };

  const handleReheat = () => {
    graphRef.current?.d3ReheatSimulation?.();
  };

  const buttonStyle: React.CSSProperties = {
    padding: "6px 10px",
    border: "1px solid #ddd",
    borderRadius: 4,
    cursor: "pointer",
    background: "white",
    fontSize: 14,
    display: "flex",
    alignItems: "center",
    justifyContent: "center",
    minWidth: 32,
  };

  return (
    <div
      style={{
        position: "absolute",
        bottom: 16,
        right: 16,
        display: "flex",
        flexDirection: "column",
        gap: 4,
        background: "white",
        borderRadius: 6,
        boxShadow: "0 2px 8px rgba(0,0,0,0.15)",
        padding: 4,
        zIndex: 50,
      }}
    >
      <button onClick={handleZoomIn} style={buttonStyle} title="Zoom in">
        +
      </button>
      <button onClick={handleZoomOut} style={buttonStyle} title="Zoom out">
        −
      </button>
      <button onClick={handleFit} style={buttonStyle} title="Fit to screen">
        ⊡
      </button>
      <button onClick={handleReheat} style={buttonStyle} title="Re-run layout simulation">
        ↻
      </button>
    </div>
  );
}
