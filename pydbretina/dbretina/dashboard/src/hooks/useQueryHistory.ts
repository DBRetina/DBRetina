import { useState, useCallback, useEffect } from "react";

export interface QueryHistoryItem {
  query: string;
  timestamp: number;
  type: "sql" | "cypher";
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

export function useQueryHistory(type: "sql" | "cypher") {
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

// Cypher query templates
// Note: Group is a reserved word in KuzuDB, so it must be backtick-escaped
export const CYPHER_TEMPLATES = [
  {
    name: "Neighbors of a Group",
    query: `MATCH (a:\`Group\`)-[r:SIMILAR_TO]->(b:\`Group\`)
WHERE a.name = "group_name"
RETURN b.name AS neighbor, r.ochiai AS similarity
ORDER BY r.ochiai DESC
LIMIT 20`,
  },
  {
    name: "Shortest Path",
    query: `MATCH (a:\`Group\`), (b:\`Group\`)
WHERE a.name = "group_a" AND b.name = "group_b"
MATCH path = shortestPath((a)-[:SIMILAR_TO*]-(b))
RETURN length(path) AS path_length`,
  },
  {
    name: "Hub Nodes (High Degree)",
    query: `MATCH (a:\`Group\`)-[r:SIMILAR_TO]-(b:\`Group\`)
WITH a.name AS group_name, COUNT(r) AS degree
WHERE degree > 10
RETURN group_name, degree
ORDER BY degree DESC
LIMIT 20`,
  },
  {
    name: "Strong Triangles",
    query: `MATCH (a:\`Group\`)-[r1:SIMILAR_TO]->(b:\`Group\`)-[r2:SIMILAR_TO]->(c:\`Group\`)-[r3:SIMILAR_TO]->(a)
WHERE r1.ochiai > 50 AND r2.ochiai > 50 AND r3.ochiai > 50
RETURN a.name, b.name, c.name,
       r1.ochiai + r2.ochiai + r3.ochiai AS total_similarity
ORDER BY total_similarity DESC
LIMIT 10`,
  },
];
