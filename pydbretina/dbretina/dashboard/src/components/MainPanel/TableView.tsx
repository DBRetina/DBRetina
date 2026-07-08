import { useState, useEffect } from "react";
import { useDashboard } from "../../state/context";
import { fetchPairs, getAPIError } from "../../api";
import { EmptyState, LoadingState } from "../Common";

export default function TableView() {
  const { state } = useDashboard();
  const [pairs, setPairs] = useState<Record<string, unknown>[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!state.info) return;
    setLoading(true);
    setError(null);

    fetchPairs(state.metric, state.cutoff, 500)
      .then((res) => {
        setPairs(res.pairs);
        setError(null);
      })
      .catch((e) => {
        setPairs([]);
        const apiError = getAPIError(e);
        setError(apiError?.detail || (e instanceof Error ? e.message : "Failed to load data"));
      })
      .finally(() => setLoading(false));
  }, [state.info, state.metric, state.cutoff]);

  if (loading) {
    return (
      <div className="table-container">
        <LoadingState message="Loading pairs..." />
      </div>
    );
  }

  if (error) {
    return (
      <div className="table-container">
        <EmptyState
          type="error"
          title="Failed to load pairs"
          message={error}
        />
      </div>
    );
  }

  if (pairs.length === 0) {
    return (
      <div className="table-container">
        <EmptyState
          type="no-data"
          title="No pairs found"
          message="No pairs meet the current filter criteria. Try lowering the cutoff value."
        />
      </div>
    );
  }

  const columns = ["group_1_name", "group_2_name", "shared_features", "containment", "ochiai", "jaccard", "csi", "dice"];

  return (
    <div className="table-container">
      <table className="pairs-table">
        <thead>
          <tr>
            {columns.map((col) => (
              <th key={col}>{col.replace(/_/g, " ")}</th>
            ))}
          </tr>
        </thead>
        <tbody>
          {pairs.map((row, i) => (
            <tr key={i}>
              {columns.map((col) => (
                <td key={col}>
                  {typeof row[col] === "number"
                    ? (row[col] as number) % 1 === 0
                      ? (row[col] as number).toLocaleString()
                      : (row[col] as number).toFixed(2)
                    : String(row[col] ?? "")}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
