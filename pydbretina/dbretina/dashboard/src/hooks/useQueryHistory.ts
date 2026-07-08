import { useState, useCallback, useEffect } from "react";

export interface QueryHistoryItem {
  query: string;
  timestamp: number;
  type: "sql";
  success: boolean;
}

const MAX_HISTORY_SIZE = 50;
const STORAGE_KEY = "dbretina_query_history";

function loadHistory(): QueryHistoryItem[] {
  try {
    const stored = localStorage.getItem(STORAGE_KEY);
    if (stored) {
      return JSON.parse(stored);
    }
  } catch {
    // Ignore parsing errors
  }
  return [];
}

function saveHistory(history: QueryHistoryItem[]) {
  try {
    localStorage.setItem(STORAGE_KEY, JSON.stringify(history.slice(0, MAX_HISTORY_SIZE)));
  } catch {
    // Ignore storage errors
  }
}

export function useQueryHistory(type: "sql") {
  const [history, setHistory] = useState<QueryHistoryItem[]>([]);

  // Load history on mount
  useEffect(() => {
    const allHistory = loadHistory();
    setHistory(allHistory.filter((h) => h.type === type));
  }, [type]);

  const addToHistory = useCallback(
    (query: string, success: boolean) => {
      const newItem: QueryHistoryItem = {
        query: query.trim(),
        timestamp: Date.now(),
        type,
        success,
      };

      setHistory((prev) => {
        // Avoid duplicates (same query)
        const filtered = prev.filter((h) => h.query !== newItem.query);
        const updated = [newItem, ...filtered].slice(0, MAX_HISTORY_SIZE);

        // Save all history
        const allHistory = loadHistory().filter((h) => h.type !== type);
        saveHistory([...updated, ...allHistory]);

        return updated;
      });
    },
    [type]
  );

  const clearHistory = useCallback(() => {
    setHistory([]);
    const allHistory = loadHistory().filter((h) => h.type !== type);
    saveHistory(allHistory);
  }, [type]);

  return { history, addToHistory, clearHistory };
}

// SQL query templates
export const SQL_TEMPLATES = [
  {
    name: "Top Similar Pairs",
    query: `SELECT group_1_name, group_2_name, ochiai, jaccard, shared_features
FROM pairs
WHERE ochiai > 50
ORDER BY ochiai DESC
LIMIT 100`,
  },
  {
    name: "Pairs by Group",
    query: `SELECT group_1_name, group_2_name, ochiai, shared_features
FROM pairs
WHERE group_1_name LIKE '%disease_name%'
   OR group_2_name LIKE '%disease_name%'
ORDER BY ochiai DESC
LIMIT 50`,
  },
  {
    name: "High Overlap Summary",
    query: `SELECT
  COUNT(*) as pair_count,
  AVG(ochiai) as avg_ochiai,
  AVG(jaccard) as avg_jaccard,
  AVG(shared_features) as avg_shared
FROM pairs
WHERE ochiai > 70`,
  },
  {
    name: "Metric Distribution",
    query: `SELECT
  CASE
    WHEN ochiai < 20 THEN '0-20'
    WHEN ochiai < 40 THEN '20-40'
    WHEN ochiai < 60 THEN '40-60'
    WHEN ochiai < 80 THEN '60-80'
    ELSE '80-100'
  END as ochiai_range,
  COUNT(*) as count
FROM pairs
GROUP BY ochiai_range
ORDER BY ochiai_range`,
  },
];
