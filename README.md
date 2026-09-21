# Swara Finder

A single-page web app (no build step) that listens to your voice or an uploaded
recording and reports the swaras (Sa Re Ga Ma Pa Dha Ni) relative to whichever
tonic (Sa) you set.

## What's in here
- `index.html` — the entire app (UI + the DSP engine from `dsp.js` inlined into it).
- `dsp.js` — kept alongside for reference/editing; `index.html` already has its
  contents baked in, so you don't need to load it separately.
- `functions/api/audio.js` — a Cloudflare Pages Function that lets you paste a
  **direct link to an audio/video file** and have the page load it, working
  around the browser's CORS restrictions.

## Deploy on Cloudflare Pages
1. Push this folder to a GitHub repo.
2. In Cloudflare Pages: **Create a project → connect the repo**.
3. Build settings: **no build command**, output directory `/` (root).
4. Deploy. The `functions/api/audio.js` file is auto-detected and deployed as
   a Pages Function — no extra config needed.

## Accuracy notes (why this should match/beat the Android app)
- **Live mic**: uses the MPM (McLeod Pitch Method) autocorrelation-in-frequency-
  domain detector in `dsp.js`, which is much more octave-stable than a naive
  autocorrelation/YIN pass — this is the same class of algorithm serious pitch
  trackers use.
- **Uploaded recordings, solo voice/instrument**: same MPM detector run frame
  by frame, then cleaned up (octave-error correction, median smoothing, gap
  bridging, tiny-blip removal) before being turned into swara segments.
- **Uploaded recordings that are a mix** (voice + tanpura/harmonium/tabla):
  a harmonic-salience engine picks out the most likely melodic line rather
  than just the loudest pitch, and — if the file is stereo — down-weights
  instruments that are hard-panned away from the centre (voice is usually
  closest to dead-centre in most recordings/mixes), plus a drone-suppression
  pass that strips out whatever pitch is *constantly* present (the tanpura
  drone) since a real melody moves and a drone doesn't.
- **Sa auto-detection**: a circular pitch-class histogram matched against a
  Hindustani/Carnatic scale-degree template, so you don't have to know your own
  tonic.

None of this is true "AI vocal isolation" (like Spotify/Adobe's ML source
separators) — it's signal-processing heuristics that work well on typical
solo-voice-forward recordings but won't cleanly separate a loud, centre-panned
harmonium from an equally loud, centre-panned voice.

## About YouTube links
The app intentionally does **not** extract audio from YouTube (or Spotify/
SoundCloud) links. Doing that means bypassing the platform's own access
controls to their streams, which their Terms of Service prohibit, and the
extraction techniques that make it "work" break constantly and are legally
gray at best. The pasted-URL box works for a **direct link to an audio/video
file** (something ending in `.mp3`, `.wav`, `.m4a`, `.mp4`, etc., hosted
anywhere). For a YouTube video: download the audio yourself through a
legitimate route (e.g. YouTube's own official download/offline feature, or
recording your own performance), then use the "drop a file" option here —
it'll get the same analysis.
