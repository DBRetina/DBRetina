/**
 * Constants for ForceGraphRenderer
 */

export const COMMUNITY_COLORS = [
  "#4361ee", "#f72585", "#4cc9f0", "#7209b7", "#f77f00",
  "#06d6a0", "#ef476f", "#118ab2", "#ffd166", "#073b4c",
  "#8338ec", "#ff6b6b", "#48bfe3", "#9b5de5", "#00f5d4",
];

export const THRESHOLDS = {
  /** Switch to 3D mode above this node count */
  USE_3D: 3000,
  /** Disable node labels above this node count */
  DISABLE_LABELS: 5000,
  /** Reduce physics simulation above this node count */
  REDUCE_PHYSICS: 10000,
  /** Fall back to Sigma for massive graphs */
  MASSIVE_GRAPH: 100000,
};

export const PHYSICS = {
  /** Alpha decay rate (higher = faster stabilization) */
  ALPHA_DECAY: 0.035,
  /** Velocity decay (higher = more friction) */
  VELOCITY_DECAY: 0.4,
  /** Warmup ticks before rendering */
  WARMUP_TICKS: 100,
  /** Cooldown ticks after interaction */
  COOLDOWN_TICKS: 200,
  /** Charge (repulsion) strength - more negative = stronger push */
  CHARGE_STRENGTH: -80,
  /** Max distance for charge interaction */
  CHARGE_DISTANCE_MAX: 500,
  /** Link distance (ideal spring length) */
  LINK_DISTANCE: 50,
  /** Link strength */
  LINK_STRENGTH: 0.3,
  /** Center gravity strength */
  CENTER_STRENGTH: 0.05,
  /** Collision radius multiplier */
  COLLIDE_RADIUS_MULT: 1.2,
};

export const COLORS = {
  PATH_HIGHLIGHT: "#f72585",
  SELECTION_HIGHLIGHT: "#4361ee",
  NODE_BORDER: "#fff",
  LINK_DEFAULT: "rgba(150,150,150,0.4)",
  LINK_PATH: "#f72585",
  TEXT: "#212529",
};
