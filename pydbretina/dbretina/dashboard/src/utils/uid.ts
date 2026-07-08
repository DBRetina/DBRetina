/**
 * Generate a UUID v4 that works in insecure contexts.
 *
 * `crypto.randomUUID()` is only available in a secure context (HTTPS or localhost).
 * `DBRetina serve` runs over plain HTTP and is typically reached via a LAN/Tailscale IP,
 * which is an insecure context, so `crypto.randomUUID` is undefined there. `crypto.getRandomValues`
 * is available in insecure contexts, so we build the UUID from it (falling back to Math.random
 * only if Web Crypto is entirely absent).
 */
export function genId(): string {
  const c: Crypto | undefined = globalThis.crypto;
  if (c && typeof c.randomUUID === "function") {
    return c.randomUUID();
  }
  const bytes = new Uint8Array(16);
  if (c && typeof c.getRandomValues === "function") {
    c.getRandomValues(bytes);
  } else {
    for (let i = 0; i < 16; i++) bytes[i] = Math.floor(Math.random() * 256);
  }
  bytes[6] = (bytes[6] & 0x0f) | 0x40; // version 4
  bytes[8] = (bytes[8] & 0x3f) | 0x80; // variant 10
  const hex = Array.from(bytes, (b) => b.toString(16).padStart(2, "0"));
  return (
    hex.slice(0, 4).join("") + "-" +
    hex.slice(4, 6).join("") + "-" +
    hex.slice(6, 8).join("") + "-" +
    hex.slice(8, 10).join("") + "-" +
    hex.slice(10, 16).join("")
  );
}
