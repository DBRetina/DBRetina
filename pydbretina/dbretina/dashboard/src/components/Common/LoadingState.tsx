interface LoadingStateProps {
  /** Message to display while loading */
  message?: string;
  /** Optional progress value (0-100) */
  progress?: number;
  /** Size of the spinner */
  size?: "small" | "medium" | "large";
  /** Whether to show as overlay on existing content */
  overlay?: boolean;
}

const spinnerSizes = {
  small: 16,
  medium: 32,
  large: 48,
};

/**
 * Loading indicator with optional message and progress bar.
 * Can be displayed inline or as an overlay.
 */
export default function LoadingState({
  message = "Loading...",
  progress,
  size = "medium",
  overlay = false,
}: LoadingStateProps) {
  const spinnerSize = spinnerSizes[size];

  const content = (
    <div
      className="loading-state"
      style={{
        display: "flex",
        flexDirection: "column",
        alignItems: "center",
        justifyContent: "center",
        gap: 12,
        padding: overlay ? 0 : 32,
      }}
    >
      {/* Spinner */}
      <div
        className="loading-spinner"
        style={{
          width: spinnerSize,
          height: spinnerSize,
          border: `${Math.max(2, spinnerSize / 12)}px solid var(--border-color)`,
          borderTopColor: "var(--accent)",
          borderRadius: "50%",
          animation: "spin 0.8s linear infinite",
        }}
      />

      {/* Message */}
      {message && (
        <div
          className="loading-message"
          style={{
            fontSize: size === "small" ? 11 : 13,
            color: "var(--text-secondary)",
          }}
        >
          {message}
        </div>
      )}

      {/* Progress bar */}
      {progress !== undefined && (
        <div
          className="loading-progress"
          style={{
            width: "100%",
            maxWidth: 200,
            height: 4,
            background: "var(--border-color)",
            borderRadius: 2,
            overflow: "hidden",
          }}
        >
          <div
            className="loading-progress-bar"
            style={{
              width: `${Math.min(100, Math.max(0, progress))}%`,
              height: "100%",
              background: "var(--accent)",
              transition: "width 0.3s ease",
            }}
          />
        </div>
      )}

      {progress !== undefined && (
        <div
          className="loading-percent"
          style={{
            fontSize: 11,
            color: "var(--text-secondary)",
          }}
        >
          {Math.round(progress)}%
        </div>
      )}
    </div>
  );

  if (overlay) {
    return (
      <div
        className="loading-overlay"
        style={{
          position: "absolute",
          inset: 0,
          display: "flex",
          alignItems: "center",
          justifyContent: "center",
          background: "rgba(255, 255, 255, 0.85)",
          zIndex: 100,
        }}
      >
        {content}
      </div>
    );
  }

  return content;
}

// Add keyframes for spinner animation
const styleSheet = document.createElement("style");
styleSheet.textContent = `
  @keyframes spin {
    to { transform: rotate(360deg); }
  }
`;
if (typeof document !== "undefined" && !document.getElementById("loading-state-styles")) {
  styleSheet.id = "loading-state-styles";
  document.head.appendChild(styleSheet);
}
