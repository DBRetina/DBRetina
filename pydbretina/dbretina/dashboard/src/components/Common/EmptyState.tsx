import { ReactNode } from "react";

interface EmptyStateProps {
  /** Title displayed prominently */
  title: string;
  /** Optional description/explanation (alias: message) */
  description?: string;
  /** Alias for description */
  message?: string;
  /** Custom icon element */
  icon?: ReactNode;
  /** Type of empty state for default icon selection */
  type?: "no-data" | "no-results" | "no-selection" | "error";
  /** Optional action button */
  action?: {
    label: string;
    onClick: () => void;
  };
  /** Size variant */
  size?: "small" | "medium" | "large";
  /** Compact mode (smaller padding, no icon) */
  compact?: boolean;
}

const defaultIcons: Record<string, ReactNode> = {
  "no-data": (
    <svg width="48" height="48" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5">
      <path d="M20 6L9 17l-5-5" />
      <circle cx="12" cy="12" r="10" strokeDasharray="4 4" />
    </svg>
  ),
  "no-results": (
    <svg width="48" height="48" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5">
      <circle cx="11" cy="11" r="8" />
      <path d="M21 21l-4.35-4.35" />
      <path d="M8 11h6" />
    </svg>
  ),
  "no-selection": (
    <svg width="48" height="48" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5">
      <circle cx="12" cy="12" r="10" />
      <path d="M12 8v4l3 3" />
    </svg>
  ),
  "error": (
    <svg width="48" height="48" viewBox="0 0 24 24" fill="none" stroke="var(--danger)" strokeWidth="1.5">
      <circle cx="12" cy="12" r="10" />
      <line x1="15" y1="9" x2="9" y2="15" />
      <line x1="9" y1="9" x2="15" y2="15" />
    </svg>
  ),
};

/**
 * Display an empty state with icon, title, description, and optional action.
 * Use when there's no data to display or no selection has been made.
 */
export default function EmptyState({
  title,
  description,
  message,
  icon,
  type = "no-data",
  action,
  size = "medium",
  compact = false,
}: EmptyStateProps) {
  // Support both description and message props
  const displayMessage = description || message;

  const sizeStyles = {
    small: { padding: "16px", iconScale: 0.75 },
    medium: { padding: "32px", iconScale: 1 },
    large: { padding: "48px", iconScale: 1.25 },
  };

  const style = compact
    ? { padding: "12px", iconScale: 0.5 }
    : sizeStyles[size];

  return (
    <div
      className="empty-state"
      style={{
        display: "flex",
        flexDirection: "column",
        alignItems: "center",
        justifyContent: "center",
        padding: style.padding,
        textAlign: "center",
        color: "var(--text-secondary)",
        height: compact ? "auto" : "100%",
      }}
    >
      {!compact && (
        <div
          className="empty-state-icon"
          style={{
            transform: `scale(${style.iconScale})`,
            marginBottom: 16,
            opacity: 0.5,
          }}
        >
          {icon || defaultIcons[type]}
        </div>
      )}

      <h3
        className="empty-state-title"
        style={{
          fontSize: compact ? 13 : (size === "small" ? 14 : 16),
          fontWeight: 600,
          marginBottom: displayMessage ? 8 : 0,
          color: "var(--text-primary)",
        }}
      >
        {title}
      </h3>

      {displayMessage && (
        <p
          className="empty-state-description"
          style={{
            fontSize: compact ? 11 : (size === "small" ? 12 : 13),
            maxWidth: 300,
            lineHeight: 1.5,
          }}
        >
          {displayMessage}
        </p>
      )}

      {action && !compact && (
        <button
          className="empty-state-action"
          onClick={action.onClick}
          style={{
            marginTop: 16,
            padding: "8px 16px",
            fontSize: 13,
            fontWeight: 500,
            color: "var(--accent)",
            background: "transparent",
            border: "1px solid var(--accent)",
            borderRadius: 4,
            cursor: "pointer",
          }}
        >
          {action.label}
        </button>
      )}
    </div>
  );
}
