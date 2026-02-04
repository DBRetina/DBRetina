import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { DashboardProvider } from "./state/context";
import App from "./App";
import "./styles/global.css";

createRoot(document.getElementById("root")!).render(
  <StrictMode>
    <DashboardProvider>
      <App />
    </DashboardProvider>
  </StrictMode>
);
