// /functions/api/audio.js
// Cloudflare Pages Function. Lets the page load a direct audio/video-file
// URL that the browser can't fetch itself because of CORS, by fetching it
// server-side and streaming the bytes back with permissive CORS headers.
//
// It deliberately does NOT extract audio from YouTube (or similar streaming
// platforms). Pulling a media stream out of YouTube means reverse-engineering
// their player and bypassing the access controls they put on their content,
// which their Terms of Service prohibit — so this endpoint returns a clear
// error for those links instead of attempting it. Everything else (a direct
// link to an .mp3/.wav/.m4a/.mp4 file hosted anywhere reachable on the public
// internet) is proxied normally.

const BLOCKED_HOSTS = [
  /(^|\.)youtube\.com$/i,
  /(^|\.)youtu\.be$/i,
  /(^|\.)youtube-nocookie\.com$/i,
  /(^|\.)music\.youtube\.com$/i,
  /(^|\.)spotify\.com$/i,
  /(^|\.)soundcloud\.com$/i, // has its own ToS-restricted stream endpoints too
];

function corsHeaders(extra) {
  return Object.assign(
    {
      "Access-Control-Allow-Origin": "*",
      "Access-Control-Allow-Methods": "GET, HEAD, OPTIONS",
      "Access-Control-Allow-Headers": "Content-Type, Range",
    },
    extra || {}
  );
}

function errorJson(status, message) {
  return new Response(JSON.stringify({ error: message }), {
    status,
    headers: corsHeaders({ "Content-Type": "application/json" }),
  });
}

// Basic SSRF guard: block anything that isn't a normal public http(s) host —
// no localhost, no bare IPs in private ranges, no link-local/metadata addresses.
function isPrivateHost(hostname) {
  const h = hostname.toLowerCase();
  if (h === "localhost" || h.endsWith(".local")) return true;
  const ipv4 = h.match(/^(\d{1,3})\.(\d{1,3})\.(\d{1,3})\.(\d{1,3})$/);
  if (ipv4) {
    const [a, b] = [parseInt(ipv4[1], 10), parseInt(ipv4[2], 10)];
    if (a === 127 || a === 10 || a === 0) return true;
    if (a === 169 && b === 254) return true; // link-local / cloud metadata
    if (a === 172 && b >= 16 && b <= 31) return true;
    if (a === 192 && b === 168) return true;
  }
  if (h === "::1" || h.startsWith("fe80:") || h.startsWith("fc") || h.startsWith("fd")) return true;
  return false;
}

export async function onRequestOptions() {
  return new Response(null, { status: 204, headers: corsHeaders() });
}

export async function onRequestGet(context) {
  const { request } = context;
  const reqUrl = new URL(request.url);
  const target = reqUrl.searchParams.get("url");

  if (!target) return errorJson(400, "Missing ?url= parameter.");

  let parsed;
  try {
    parsed = new URL(target);
  } catch (e) {
    return errorJson(400, "That is not a valid URL.");
  }

  if (parsed.protocol !== "https:" && parsed.protocol !== "http:") {
    return errorJson(400, "Only http/https links are supported.");
  }
  if (isPrivateHost(parsed.hostname)) {
    return errorJson(400, "That host can't be fetched.");
  }
  if (BLOCKED_HOSTS.some((re) => re.test(parsed.hostname))) {
    return errorJson(
      422,
      "This site can't pull audio out of YouTube/Spotify/SoundCloud links directly — extracting a stream from those platforms means bypassing their access controls, which isn't something this tool does. Download the audio yourself using the platform's own export/download feature (or record it) and upload the file instead, or paste a direct link to an audio/video file hosted elsewhere."
    );
  }

  let upstream;
  try {
    upstream = await fetch(parsed.toString(), {
      method: "GET",
      redirect: "follow",
      headers: {
        // A plain UA/Accept helps a few CDNs that reject requests with no
        // browser-like headers; this is a normal file fetch, not spoofing.
        "User-Agent": "Mozilla/5.0 (compatible; SwaraFinderAudioProxy/1.0)",
        Accept: "audio/*,video/*;q=0.9,*/*;q=0.5",
      },
    });
  } catch (e) {
    return errorJson(502, "Could not reach that URL.");
  }

  if (!upstream.ok) {
    return errorJson(upstream.status, "The link responded with status " + upstream.status + ".");
  }

  const ct = (upstream.headers.get("content-type") || "").toLowerCase();
  if (ct.includes("text/html")) {
    return errorJson(
      422,
      "That link points to a web page, not an audio/video file. Find the page's direct media file link (often ending in .mp3, .wav, .m4a, or .mp4) and paste that instead."
    );
  }

  const headers = corsHeaders({
    "Content-Type": ct || "application/octet-stream",
    "Cache-Control": "public, max-age=3600",
  });
  const len = upstream.headers.get("content-length");
  if (len) headers["Content-Length"] = len;

  return new Response(upstream.body, { status: 200, headers });
}
