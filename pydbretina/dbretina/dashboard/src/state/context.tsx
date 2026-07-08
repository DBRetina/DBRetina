import { createContext, useContext, useReducer, ReactNode, Dispatch } from "react";
import { DashboardState, DashboardAction } from "./types";
import { dashboardReducer, initialState } from "./reducer";

const DashboardContext = createContext<{
  state: DashboardState;
  dispatch: Dispatch<DashboardAction>;
} | null>(null);

export function DashboardProvider({ children }: { children: ReactNode }) {
  const [state, dispatch] = useReducer(dashboardReducer, initialState);
  return (
    <DashboardContext.Provider value={{ state, dispatch }}>
      {children}
    </DashboardContext.Provider>
  );
}

export function useDashboard() {
  const ctx = useContext(DashboardContext);
  if (!ctx) throw new Error("useDashboard must be used within DashboardProvider");
  return ctx;
}
