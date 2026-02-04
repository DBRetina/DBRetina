/**
 * Types for react-force-graph integration
 */

export interface ForceNode {
  id: string;
  label: string;
  degree: number;
  community: number;
  pagerank: number;
  // Visual properties
  color: string;
  size: number;
  // Computed state
  isHighlighted: boolean;
  isSelected: boolean;
  isOnPath: boolean;
  // Force simulation positions (set by d3-force)
  x?: number;
  y?: number;
  z?: number;
  // Fixed positions (to pin nodes)
  fx?: number;
  fy?: number;
  fz?: number;
  // Velocity (internal use)
  vx?: number;
  vy?: number;
  vz?: number;
}

export interface ForceLink {
  source: string | ForceNode;
  target: string | ForceNode;
  weight: number;
  shared_features: number;
  // Visual properties
  color: string;
  width: number;
  isOnPath: boolean;
}

export interface ForceGraphData {
  nodes: ForceNode[];
  links: ForceLink[];
}
