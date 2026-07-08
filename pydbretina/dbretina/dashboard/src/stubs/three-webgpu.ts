// Stub for the `three/webgpu` entry point.
//
// `three-render-objects` does `import { WebGPURenderer } from "three/webgpu"` for an
// opt-in WebGPU renderer that DBRetina never enables (it renders with WebGLRenderer;
// `useWebGPU` defaults to false and the 3D path is gated behind a 3000-node threshold).
//
// three's real `three/webgpu` build crashes at module init on WebGPU-capable browsers
// (a circular dependency leaves `NodeShaderStage` undefined -> "Cannot read properties
// of undefined (reading 'VERTEX')"), which white-screens the whole dashboard on devices
// that expose `navigator.gpu`. Aliasing the specifier to this stub keeps that module out
// of the bundle entirely while leaving WebGLRenderer untouched.
export class WebGPURenderer {
  constructor() {
    throw new Error("three/webgpu is stubbed out in the DBRetina dashboard (WebGL is used).");
  }
}
