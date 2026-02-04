import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  base: "/",
  build: {
    outDir: "../dashboard_dist",
    emptyOutDir: true,
    commonjsOptions: {
      include: [/lodash-es/, /node_modules/],
    },
  },
  optimizeDeps: {
    include: [
      "lodash-es",
      "force-graph",
      "react-force-graph-2d",
      "react-force-graph-3d",
      "three",
    ],
  },
  resolve: {
    dedupe: ["lodash-es", "three"],
  },
  server: {
    port: 5173,
    proxy: {
      "/api": "http://localhost:8080",
    },
  },
});
