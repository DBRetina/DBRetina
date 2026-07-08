import { useState, useEffect } from "react";
import { useDashboard } from "../../state/context";
import {
  BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer,
  ScatterChart, Scatter,
} from "recharts";
import { fetchTopPairs, fetchMetricSummary } from "../../api";

export default function StatsView() {
  const { state } = useDashboard();
  const [topPairs, setTopPairs] = useState<Record<string, unknown>[]>([]);
  const [summary, setSummary] = useState<Record<string, number> | null>(null);

  useEffect(() => {
    if (!state.info) return;
    fetchTopPairs(state.metric, 20).then((res) => setTopPairs(res.pairs)).catch(() => {});
    fetchMetricSummary(state.metric).then(setSummary).catch(() => {});
  }, [state.info, state.metric]);

  // Degree distribution from graph data
  const degreeDist = (() => {
    if (!state.graphData) return [];
    const counts: Record<number, number> = {};
    for (const node of state.graphData.nodes) {
      counts[node.degree] = (counts[node.degree] || 0) + 1;
    }
    return Object.entries(counts)
      .map(([deg, count]) => ({ degree: Number(deg), count }))
      .sort((a, b) => a.degree - b.degree);
  })();

  const topPairsData = topPairs.slice(0, 20).map((p: any, i: number) => ({
    name: `${String(p.group_1_name || p.group_1_id).slice(0, 25)} - ${String(p.group_2_name || p.group_2_id).slice(0, 25)}`,
    value: Number(p[state.metric]) || 0,
    idx: i,
  })).reverse();

  return (
    <div className="stats-grid">
      <div className="stats-card">
        <h3>Metric Summary: {state.metric}</h3>
        {summary && (
          <div style={{ fontSize: 13 }}>
            {Object.entries(summary).map(([k, v]) => (
              <div key={k} className="detail-row">
                <span className="detail-label">{k}</span>
                <span className="detail-value">{typeof v === "number" ? v.toFixed(2) : String(v)}</span>
              </div>
            ))}
          </div>
        )}
      </div>

      <div className="stats-card">
        <h3>Degree Distribution</h3>
        <ResponsiveContainer width="100%" height={220}>
          <ScatterChart>
            <CartesianGrid strokeDasharray="3 3" />
            <XAxis dataKey="degree" name="Degree" type="number" scale="log" domain={["auto", "auto"]} fontSize={10} />
            <YAxis dataKey="count" name="Count" type="number" scale="log" domain={["auto", "auto"]} fontSize={10} />
            <Tooltip />
            <Scatter data={degreeDist} fill="#4361ee" />
          </ScatterChart>
        </ResponsiveContainer>
      </div>

      <div className="stats-card" style={{ gridColumn: "1 / -1" }}>
        <h3>Top 20 Pairs by {state.metric}</h3>
        <ResponsiveContainer width="100%" height={Math.max(250, topPairsData.length * 22)}>
          <BarChart data={topPairsData} layout="vertical" margin={{ left: 200 }}>
            <CartesianGrid strokeDasharray="3 3" />
            <XAxis type="number" fontSize={10} />
            <YAxis type="category" dataKey="name" width={200} fontSize={9} />
            <Tooltip />
            <Bar dataKey="value" fill="#4361ee" />
          </BarChart>
        </ResponsiveContainer>
      </div>
    </div>
  );
}
