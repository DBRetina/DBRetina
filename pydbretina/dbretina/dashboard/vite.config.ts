import { fileURLToPath, URL } from "node:url";
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
    alias: {
      // three/webgpu crashes at init on WebGPU browsers and is never used here; stub it out.
      "three/webgpu": fileURLToPath(new URL("./src/stubs/three-webgpu.ts", import.meta.url)),
    },
  },
  server: {
    port: 5173,
    proxy: {
      "/api": "http://localhost:8080",
    },
  },
});
