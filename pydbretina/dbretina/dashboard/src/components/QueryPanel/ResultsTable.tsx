import { QueryResult } from "../../state/types";

interface Props {
  result: QueryResult;
}

export default function ResultsTable({ result }: Props) {
  if (result.row_count === 0) {
    return (
      <div style={{ padding: "8px 12px", fontSize: 12, color: "var(--text-secondary)" }}>
        No results
      </div>
    );
  }

  return (
    <div className="query-results">
      <table>
        <thead>
          <tr>
            {result.columns.map((col) => (
              <th key={col}>{col}</th>
            ))}
          </tr>
        </thead>
        <tbody>
          {result.rows.map((row, i) => (
            <tr key={i}>
              {result.columns.map((col) => (
                <td key={col}>
                  {typeof row[col] === "number"
                    ? (row[col] as number) % 1 === 0
                      ? (row[col] as number).toLocaleString()
                      : (row[col] as number).toFixed(4)
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
