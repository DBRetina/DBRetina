import { useState, useRef, ReactNode } from "react";

interface CollapsibleSectionProps {
  title: string;
  children: ReactNode;
  defaultExpanded?: boolean;
  loading?: boolean;
  error?: string | null;
  onExpand?: () => void;
  onRetry?: () => void;
  badge?: string | number;
}

export default function CollapsibleSection({
  title,
  children,
  defaultExpanded = false,
  loading = false,
  error = null,
  onExpand,
  onRetry,
  badge,
}: CollapsibleSectionProps) {
  const [expanded, setExpanded] = useState(defaultExpanded);
  const hasExpanded = useRef(defaultExpanded);

  function handleToggle() {
    const willExpand = !expanded;
    setExpanded(willExpand);
    if (willExpand && !hasExpanded.current) {
      hasExpanded.current = true;
      onExpand?.();
    }
  }

  return (
    <div className="detail-section">
      <div
        onClick={handleToggle}
        style={{
          display: "flex",
          alignItems: "center",
          justifyContent: "space-between",
          cursor: "pointer",
          userSelect: "none",
        }}
      >
        <h4 style={{ margin: 0, display: "flex", alignItems: "center", gap: 6 }}>
          {title}
          {badge !== undefined && (
            <span
              style={{
                fontSize: 10,
                background: "var(--accent)",
                color: "white",
                padding: "1px 6px",
                borderRadius: 8,
                fontWeight: 500,
              }}
            >
              {badge}
            </span>
          )}
        </h4>
        <span style={{ fontSize: 12, color: "var(--text-secondary)" }}>
          {expanded ? "[-]" : "[+]"}
        </span>
      </div>

      {expanded && (
        <div style={{ marginTop: 8 }}>
          {loading && (
            <div style={{ fontSize: 11, color: "var(--text-secondary)", padding: "4px 0" }}>
              Loading...
            </div>
          )}
          {error && (
            <div style={{ fontSize: 11, color: "#e74c3c", padding: "4px 0" }}>
              {error}
              {onRetry && (
                <button
                  onClick={(e) => {
                    e.stopPropagation();
                    onRetry();
                  }}
                  style={{
                    marginLeft: 8,
                    padding: "2px 8px",
                    fontSize: 10,
                    border: "1px solid var(--border-color)",
                    borderRadius: 3,
                    background: "white",
                    cursor: "pointer",
                  }}
                >
                  Retry
                </button>
              )}
            </div>
          )}
          {!loading && !error && children}
        </div>
      )}
    </div>
  );
}
