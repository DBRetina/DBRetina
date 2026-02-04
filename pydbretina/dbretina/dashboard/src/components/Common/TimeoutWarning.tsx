import { useEffect, useState } from "react";

interface TimeoutWarningProps {
  /** When the operation started (timestamp in ms) */
  startTime: number;
  /** Threshold in seconds to show warning */
  warningThreshold?: number;
  /** Message to display */
  message?: string;
  /** Callback when user cancels */
  onCancel?: () => void;
}

/**
 * Warning component that appears when an operation is taking too long.
 * Shows elapsed time and offers a cancel option.
 */
export default function TimeoutWarning({
  startTime,
  warningThreshold = 5,
  message = "This operation is taking longer than expected",
  onCancel,
}: TimeoutWarningProps) {
  const [elapsed, setElapsed] = useState(0);
  const [visible, setVisible] = useState(false);

  useEffect(() => {
    const interval = setInterval(() => {
      const seconds = Math.floor((Date.now() - startTime) / 1000);
      setElapsed(seconds);
      if (seconds >= warningThreshold && !visible) {
        setVisible(true);
      }
    }, 1000);

    return () => clearInterval(interval);
  }, [startTime, warningThreshold, visible]);

  if (!visible) {
    return null;
  }

  const formatTime = (seconds: number): string => {
    if (seconds < 60) {
      return `${seconds}s`;
    }
    const mins = Math.floor(seconds / 60);
    const secs = seconds % 60;
    return `${mins}m ${secs}s`;
  };

  return (
    <div
      className="timeout-warning"
      style={{
        display: "flex",
        alignItems: "center",
        gap: 12,
        padding: "10px 14px",
        background: "var(--warning)",
        color: "#000",
        borderRadius: 6,
        fontSize: 13,
        animation: "fadeIn 0.3s ease",
      }}
    >
      {/* Warning icon */}
      <svg
        width="18"
        height="18"
        viewBox="0 0 24 24"
        fill="none"
        stroke="currentColor"
        strokeWidth="2"
      >
        <path d="M10.29 3.86L1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z" />
        <line x1="12" y1="9" x2="12" y2="13" />
        <line x1="12" y1="17" x2="12.01" y2="17" />
      </svg>

      <div style={{ flex: 1 }}>
        <div style={{ fontWeight: 500 }}>{message}</div>
        <div style={{ fontSize: 11, opacity: 0.8, marginTop: 2 }}>
          Elapsed: {formatTime(elapsed)}
        </div>
      </div>

      {onCancel && (
        <button
          onClick={onCancel}
          style={{
            padding: "4px 10px",
            fontSize: 12,
            fontWeight: 500,
            background: "rgba(0,0,0,0.1)",
            border: "none",
            borderRadius: 4,
            cursor: "pointer",
          }}
        >
          Cancel
        </button>
      )}
    </div>
  );
}

// Add fadeIn animation
const styleSheet = document.createElement("style");
styleSheet.textContent = `
  @keyframes fadeIn {
    from { opacity: 0; transform: translateY(-10px); }
    to { opacity: 1; transform: translateY(0); }
  }
`;
if (typeof document !== "undefined" && !document.getElementById("timeout-warning-styles")) {
  styleSheet.id = "timeout-warning-styles";
  document.head.appendChild(styleSheet);
}
