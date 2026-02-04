import { useMemo } from "react";
import { useDashboard } from "../../state/context";

const COMMUNITY_COLORS = [
  "#4361ee", "#f72585", "#4cc9f0", "#7209b7", "#f77f00",
  "#06d6a0", "#ef476f", "#118ab2", "#ffd166", "#073b4c",
  "#8338ec", "#ff6b6b", "#48bfe3", "#9b5de5", "#00f5d4",
];

export default function CommunityPanel() {
  const { state, dispatch } = useDashboard();

  const communities = useMemo(() => {
    if (!state.graphData) return [];
    const counts: Record<number, number> = {};
    for (const node of state.graphData.nodes) {
      counts[node.community] = (counts[node.community] || 0) + 1;
    }
    return Object.entries(counts)
      .sort((a, b) => b[1] - a[1])
      .map(([cid, count]) => ({ id: Number(cid), count }));
  }, [state.graphData]);

  function handleClick(communityId: number) {
    if (!state.graphData) return;
    const ids = new Set(
      state.graphData.nodes
        .filter((n) => n.community === communityId)
        .map((n) => n.id)
    );
    dispatch({ type: "SET_HIGHLIGHT_NODES", nodes: ids });
  }

  if (communities.length === 0) return null;

  return (
    <div className="sidebar-section">
      <label>Communities ({communities.length})</label>
      <div className="community-list">
        {communities.slice(0, 20).map(({ id, count }) => (
          <div key={id} className="community-item" onClick={() => handleClick(id)}>
            <span
              className="community-dot"
              style={{ background: COMMUNITY_COLORS[id % COMMUNITY_COLORS.length] }}
            />
            <span>Community {id}</span>
            <span style={{ marginLeft: "auto", color: "var(--text-secondary)" }}>
              {count}
            </span>
          </div>
        ))}
      </div>
    </div>
  );
}
