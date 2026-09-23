/* ==================== DSP CORE ==================== */
var DSP = (function(){
'use strict';
var PI2 = Math.PI * 2;

function nextPow2(n){ var p = 1; while (p < n) p <<= 1; return p; }
function hz2midi(f){ return 69 + 12 * Math.log2(f / 440); }
function midi2hz(m){ return 440 * Math.pow(2, (m - 69) / 12); }
function median(a){
  var n = a.length; if (!n) return NaN;
  var s = Array.prototype.slice.call(a).sort(function(x, y){ return x - y; });
  return n % 2 ? s[(n - 1) >> 1] : 0.5 * (s[n / 2 - 1] + s[n / 2]);
}
function now(){ return (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now(); }

/* ---------- FFT (radix-2, in place) ---------- */
function FFT(n){
  this.n = n;
  var bits = Math.round(Math.log2(n)), i;
  this.cos = new Float64Array(n >> 1);
  this.sin = new Float64Array(n >> 1);
  for (i = 0; i < (n >> 1); i++){ this.cos[i] = Math.cos(PI2 * i / n); this.sin[i] = Math.sin(PI2 * i / n); }
  this.rev = new Uint32Array(n);
  for (i = 0; i < n; i++) this.rev[i] = (this.rev[i >> 1] >> 1) | ((i & 1) << (bits - 1));
}
FFT.prototype.forward = function(re, im){
  var n = this.n, rev = this.rev, cos = this.cos, sin = this.sin, i, j, k, t;
  for (i = 0; i < n; i++){
    j = rev[i];
    if (j > i){ t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
  }
  for (var size = 2; size <= n; size <<= 1){
    var half = size >> 1, step = n / size;
    for (i = 0; i < n; i += size){
      for (j = i, k = 0; j < i + half; j++, k += step){
        var l = j + half;
        var tr = re[l] * cos[k] + im[l] * sin[k];
        var ti = -re[l] * sin[k] + im[l] * cos[k];
        re[l] = re[j] - tr; im[l] = im[j] - ti;
        re[j] += tr; im[j] += ti;
      }
    }
  }
};
FFT.prototype.inverse = function(re, im){
  this.forward(im, re);
  var s = 1 / this.n;
  for (var i = 0; i < this.n; i++){ re[i] *= s; im[i] *= s; }
};

function hann(N){
  var w = new Float32Array(N);
  for (var i = 0; i < N; i++) w[i] = 0.5 - 0.5 * Math.cos(PI2 * i / (N - 1));
  return w;
}

/* ---------- windowed-sinc decimator ---------- */
function decimate(x, factor){
  if (factor <= 1) return x;
  var M = 8 * factor, taps = 2 * M + 1, h = new Float64Array(taps), fc = 0.45 / factor, sum = 0, i, k;
  for (i = 0; i < taps; i++){
    var d = i - M;
    var s = d === 0 ? 2 * fc : Math.sin(PI2 * fc * d) / (Math.PI * d);
    h[i] = s * (0.54 - 0.46 * Math.cos(PI2 * i / (taps - 1)));
    sum += h[i];
  }
  for (i = 0; i < taps; i++) h[i] /= sum;
  var n = Math.floor(x.length / factor), out = new Float32Array(n), len = x.length;
  for (var j = 0; j < n; j++){
    var c = j * factor, acc = 0;
    if (c - M >= 0 && c + M < len){
      for (k = 0; k < taps; k++) acc += h[k] * x[c - M + k];
    } else {
      for (k = 0; k < taps; k++){ var idx = c - M + k; if (idx >= 0 && idx < len) acc += h[k] * x[idx]; }
    }
    out[j] = acc;
  }
  return out;
}

/* ---------- MPM / NSDF pitch detector (FFT based, O(N log N)) ---------- */
function MPM(sr, W, fmin, fmax, kThresh){
  var fftN = nextPow2(2 * W), fft = new FFT(fftN);
  var re = new Float64Array(fftN), im = new Float64Array(fftN);
  var xs = new Float64Array(W), nsdf = new Float32Array(W);
  var tauMin = Math.max(2, Math.floor(sr / fmax));
  var tauMax = Math.min(W - 3, Math.ceil(sr / fmin));
  kThresh = kThresh || 0.9;
  this.W = W;
  this.detect = function(x, off){
    var i, mean = 0, e = 0;
    for (i = 0; i < W; i++) mean += x[off + i];
    mean /= W;
    re.fill(0); im.fill(0);
    for (i = 0; i < W; i++){ var v = x[off + i] - mean; xs[i] = v; re[i] = v; e += v * v; }
    var rms = Math.sqrt(e / W);
    if (e < 1e-12) return { f0: 0, clarity: 0, rms: rms };
    fft.forward(re, im);
    for (i = 0; i < fftN; i++){ re[i] = re[i] * re[i] + im[i] * im[i]; im[i] = 0; }
    fft.inverse(re, im);
    var m = 2 * e;
    nsdf[0] = 1;
    for (var tau = 1; tau <= tauMax + 1 && tau < W; tau++){
      m -= xs[tau - 1] * xs[tau - 1] + xs[W - tau] * xs[W - tau];
      nsdf[tau] = m > 1e-12 ? 2 * re[tau] / m : 0;
    }
    // key maxima: highest point of each positive lobe (after the zero-lag lobe)
    var t = 1, nmax = -1, keys = [];
    while (t < tauMax && nsdf[t] > 0) t++;
    while (t < tauMax){
      while (t < tauMax && nsdf[t] <= 0) t++;
      if (t >= tauMax) break;
      var best = t;
      while (t < tauMax && nsdf[t] > 0){ if (nsdf[t] > nsdf[best]) best = t; t++; }
      if (best >= tauMin){
        keys.push(best);
        if (nsdf[best] > nmax) nmax = nsdf[best];
      }
    }
    if (!keys.length || nmax < 0.3) return { f0: 0, clarity: nmax > 0 ? nmax : 0, rms: rms };
    var thr = kThresh * nmax, pick = keys[0];
    for (i = 0; i < keys.length; i++){ if (nsdf[keys[i]] >= thr){ pick = keys[i]; break; } }
    // Bass-safe fundamental check: on low, harmonic-rich notes a strong
    // harmonic can create a spurious short-lag (upper-octave) peak that
    // clears the threshold above before the true, longer-lag fundamental
    // peak is reached (short windows only fit a couple of true periods, so
    // the fundamental's own peak can look weaker than a harmonic's).
    // If a comparably strong lobe sits at a clean integer multiple of the
    // chosen lag, it is almost certainly the real fundamental — the short
    // lag is then just that harmonic re-matching itself — so walk out to
    // the longest such well-supported multiple instead.
    for (i = 0; i < keys.length; i++){
      var kk = keys[i], ratio = kk / pick, mult = Math.round(ratio);
      if (mult >= 2 && mult <= 4 && Math.abs(ratio - mult) < 0.08 && nsdf[kk] >= 0.92 * nsdf[pick]){
        pick = kk;
      }
    }
    var a = nsdf[pick - 1], b = nsdf[pick], c = nsdf[pick + 1];
    var den = a - 2 * b + c, shift = den !== 0 ? 0.5 * (a - c) / den : 0;
    if (shift > 1 || shift < -1) shift = 0;
    var tauF = pick + shift, clarity = b - 0.25 * (a - c) * shift;
    return { f0: sr / tauF, clarity: clarity, rms: rms };
  };
}

/* ---------- Harmonic-salience engine (melody in a mixture) ----------
   Windowed FFT -> spectral peaks -> harmonic summation onto a log-frequency grid.
   If a right channel is supplied, bins are weighted by "centre-panned-ness"
   (in-phase, equal level in L and R) so hard-panned instruments count for less. */
function SalienceEngine(sr, o){
  o = o || {};
  var N = o.N || 2048, fftN = o.fftN || 4096;
  var fmin = o.fmin || 75, fmax = o.fmax || 1100, binCents = o.binCents || 20;
  var H = o.H || 10, alpha = o.alpha || 0.8, beta = o.beta || 1.0;
  var maxPeaks = o.maxPeaks || 32;
  var centerPow = o.centerPow === undefined ? 2 : o.centerPow;
  var nBins = Math.floor(1200 * Math.log2(fmax / fmin) / binCents) + 1;
  var win = hann(N), fft = new FFT(fftN);
  var re = new Float64Array(fftN), im = new Float64Array(fftN);
  var half = fftN >> 1, mag = new Float32Array(half + 1), binHz = sr / fftN;
  var kLo = Math.max(2, Math.floor(50 / binHz)), kHi = Math.min(half - 2, Math.ceil(5000 / binHz));
  var candF = new Float32Array(400), candA = new Float32Array(400);
  var hw = new Float32Array(H + 1); for (var h0 = 1; h0 <= H; h0++) hw[h0] = Math.pow(alpha, h0 - 1);
  var spread = Math.ceil(50 / binCents) + 1;
  var cwTab = new Float32Array(51);
  for (var ci = 0; ci <= 50; ci++){ var cc = Math.cos(ci / 50 * Math.PI / 2); cwTab[ci] = cc * cc; }
  var self = this;
  self.nBins = nBins; self.N = N; self.binCents = binCents; self.fmin = fmin; self.fmax = fmax;
  self.maxPeaks = maxPeaks; self.H = H; self.alpha = alpha;
  self.sal = new Float32Array(nBins);
  self.peakF = new Float32Array(maxPeaks); self.peakA = new Float32Array(maxPeaks);
  self.nPeaks = 0; self.maxA = 0; self.mass = 0; self.salMax = 0;

  self.process = function(L, R, off){
    var i, k, len = L.length;
    for (i = 0; i < N; i++){
      var p = off + i, w = win[i];
      re[i] = (p < len ? L[p] : 0) * w;
      im[i] = R ? (p < len ? R[p] : 0) * w : 0;
    }
    for (; i < fftN; i++){ re[i] = 0; im[i] = 0; }
    fft.forward(re, im);
    if (R){
      for (k = 0; k <= half; k++){
        var kk = (fftN - k) & (fftN - 1);
        var lr = (re[k] + re[kk]) * 0.5, li = (im[k] - im[kk]) * 0.5;
        var rr = (im[k] + im[kk]) * 0.5, ri = -(re[k] - re[kk]) * 0.5;
        var den = lr * lr + li * li + rr * rr + ri * ri + 1e-20;
        var psi = 2 * (lr * rr + li * ri) / den;
        var mr = lr + rr, mi = li + ri, m = 0.5 * Math.sqrt(mr * mr + mi * mi);
        mag[k] = psi > 0 ? m * (centerPow === 2 ? psi * psi : Math.pow(psi, centerPow)) : 0;
      }
    } else {
      for (k = 0; k <= half; k++) mag[k] = Math.sqrt(re[k] * re[k] + im[k] * im[k]);
    }
    var mx = 0;
    for (k = kLo; k <= kHi; k++) if (mag[k] > mx) mx = mag[k];
    self.maxA = mx;
    var sal = self.sal; sal.fill(0);
    if (mx <= 1e-9){ self.nPeaks = 0; self.mass = 0; self.salMax = 0; return; }
    var thr = mx * 0.01, nc = 0;
    for (k = kLo; k <= kHi && nc < 400; k++){
      var b = mag[k];
      if (b > thr && b > mag[k - 1] && b >= mag[k + 1]){
        var a1 = mag[k - 1], c1 = mag[k + 1], dn = a1 - 2 * b + c1;
        var d = dn !== 0 ? 0.5 * (a1 - c1) / dn : 0;
        if (d > 0.5 || d < -0.5) d = 0;
        candF[nc] = (k + d) * binHz; candA[nc] = b - 0.25 * (a1 - c1) * d; nc++;
      }
    }
    var np = Math.min(nc, maxPeaks), j, t;
    for (j = 0; j < np; j++){
      var bi = j;
      for (i = j + 1; i < nc; i++) if (candA[i] > candA[bi]) bi = i;
      if (bi !== j){ t = candA[j]; candA[j] = candA[bi]; candA[bi] = t; t = candF[j]; candF[j] = candF[bi]; candF[bi] = t; }
    }
    var mass = 0, smax = 0;
    for (j = 0; j < np; j++){
      var pa = beta === 1 ? candA[j] / mx : Math.pow(candA[j] / mx, beta), pf = candF[j];
      self.peakF[j] = pf; self.peakA[j] = candA[j];
      mass += pa;
      for (var h = 1; h <= H; h++){
        var f0 = pf / h;
        if (f0 < fmin * 0.97) break;
        if (f0 > fmax * 1.03) continue;
        var pos = 1200 * Math.log2(f0 / fmin) / binCents, c = Math.round(pos), wgt = pa * hw[h];
        for (var dd = -spread; dd <= spread; dd++){
          var bb = c + dd;
          if (bb < 0 || bb >= nBins) continue;
          var dist = Math.abs(bb - pos) * binCents;
          if (dist >= 50) continue;
          sal[bb] += wgt * cwTab[(dist + 0.5) | 0];
        }
      }
    }
    for (i = 0; i < nBins; i++) if (sal[i] > smax) smax = sal[i];
    self.nPeaks = np; self.mass = mass; self.salMax = smax;
  };

  /* refine an f0 guess with the actual spectral peaks of the current frame */
  self.refine = function(f0, pf, pa, np){
    var cur = f0;
    for (var it = 0; it < 2; it++){
      var num = 0, den = 0;
      for (var h = 1; h <= 8; h++){
        var target = cur * h, bestJ = -1, bestD = 35;
        for (var j = 0; j < np; j++){
          var dc = Math.abs(1200 * Math.log2(pf[j] / target));
          if (dc < bestD){ bestD = dc; bestJ = j; }
        }
        if (bestJ >= 0){ var wt = pa[bestJ] * hw[h]; num += wt * (pf[bestJ] / h); den += wt; }
      }
      if (den > 0) cur = num / den; else break;
    }
    return cur;
  };
}

/* voice-range prior over the salience grid */
function rangePrior(nBins, fmin, binCents, loHz, hiHz){
  var pr = new Float32Array(nBins);
  for (var b = 0; b < nBins; b++){
    var f = fmin * Math.pow(2, b * binCents / 1200), s = 1;
    if (f < loHz) s = Math.max(0.05, 1 - (1200 * Math.log2(loHz / f)) / 600);
    else if (f > hiHz) s = Math.max(0.05, 1 - (1200 * Math.log2(f / hiHz)) / 600);
    pr[b] = s;
  }
  return pr;
}

/* ---------- MELODY EXTRACTION for a mix (generator so the UI can show progress) ---------- */
function* melodyGen(L, R, sr, o){
  o = o || {};
  var N = 2048, hop = o.hop || 384;
  var loHz = o.loHz || 80, hiHz = o.hiHz || 1000;
  var eng = new SalienceEngine(sr, { N: N, fftN: 4096, fmin: 70, fmax: 1150, binCents: 20, centerPow: R ? (o.centerPow === undefined ? 2 : o.centerPow) : 0 });
  var T = Math.max(1, Math.floor((L.length - N) / hop) + 1), nBins = eng.nBins, MP = eng.maxPeaks;
  var raw = new Float32Array(T * nBins), obs = new Uint8Array(T * nBins);
  var salMax = new Float32Array(T), mass = new Float32Array(T), maxA = new Float32Array(T);
  var PF = new Float32Array(T * MP), PA = new Float32Array(T * MP), PN = new Uint8Array(T);
  var prior = rangePrior(nBins, eng.fmin, eng.binCents, loHz, hiHz);
  var t, b, last = now();
  for (t = 0; t < T; t++){
    eng.process(L, R, t * hop);
    mass[t] = eng.mass; maxA[t] = eng.maxA; PN[t] = eng.nPeaks;
    raw.set(eng.sal, t * nBins);
    for (var q = 0; q < eng.nPeaks; q++){ PF[t * MP + q] = eng.peakF[q]; PA[t * MP + q] = eng.peakA[q]; }
    if ((t & 31) === 0 && now() - last > 40){ last = now(); yield 0.6 * t / T; }
  }
  // ---- drone suppression: remove what is steady in each pitch bin over +-3 s (tanpura / shruti box) ----
  var gamma = o.drone === undefined ? 0.75 : o.drone;
  if (gamma > 0 && T > 8){
    var w = Math.max(4, Math.round(3 / (hop / sr))), pre = new Float64Array(T + 1);
    for (b = 0; b < nBins; b++){
      for (t = 0; t < T; t++) pre[t + 1] = pre[t] + raw[t * nBins + b];
      for (t = 0; t < T; t++){
        var a0 = t - w < 0 ? 0 : t - w, a1 = t + w + 1 > T ? T : t + w + 1;
        var mean = (pre[a1] - pre[a0]) / (a1 - a0), v = raw[t * nBins + b] - gamma * mean;
        raw[t * nBins + b] = v > 0 ? v : 0;
      }
      if ((b & 15) === 0 && now() - last > 40){ last = now(); yield 0.6 + 0.05 * b / nBins; }
    }
  }
  for (t = 0; t < T; t++){
    var base = t * nBins, sm = 0;
    for (b = 0; b < nBins; b++){ var vv = raw[base + b] * prior[b]; raw[base + b] = vv; if (vv > sm) sm = vv; }
    salMax[t] = sm;
    var inv = sm > 0 ? 255 / sm : 0;
    for (b = 0; b < nBins; b++) obs[base + b] = Math.round(raw[base + b] * inv);
  }
  // ---- Viterbi with linear jump penalty (two-pass distance transform) ----
  var tau = o.jumpPenalty || 0.006;
  var back = new Int16Array(T * nBins);
  var cur = new Float32Array(nBins), g = new Float32Array(nBins), s = new Int16Array(nBins);
  for (b = 0; b < nBins; b++) cur[b] = obs[b] / 255;
  for (t = 1; t < T; t++){
    for (b = 0; b < nBins; b++){ g[b] = cur[b]; s[b] = b; }
    for (b = 1; b < nBins; b++) if (g[b - 1] - tau > g[b]){ g[b] = g[b - 1] - tau; s[b] = s[b - 1]; }
    for (b = nBins - 2; b >= 0; b--) if (g[b + 1] - tau > g[b]){ g[b] = g[b + 1] - tau; s[b] = s[b + 1]; }
    var bt = t * nBins, mxv = -1e9;
    for (b = 0; b < nBins; b++){ back[bt + b] = s[b]; cur[b] = g[b] + obs[bt + b] / 255; if (cur[b] > mxv) mxv = cur[b]; }
    for (b = 0; b < nBins; b++) cur[b] -= mxv;
    if ((t & 63) === 0 && now() - last > 40){ last = now(); yield 0.6 + 0.25 * t / T; }
  }
  var path = new Int16Array(T), bestB = 0;
  for (b = 1; b < nBins; b++) if (cur[b] > cur[bestB]) bestB = b;
  path[T - 1] = bestB;
  for (t = T - 1; t > 0; t--) path[t - 1] = back[t * nBins + path[t]];
  yield 0.9;

  // ---- refine, voicing ----
  var midi = new Float32Array(T), fracArr = new Float32Array(T), sp = new Float32Array(T), pitchAll = new Float32Array(T);
  var sortedMax = Array.prototype.slice.call(maxA).sort(function(a, c){ return a - c; });
  var p90 = sortedMax[Math.floor(0.9 * (T - 1))] || 1e-9;
  var fracThr = o.voicing === undefined ? 0.28 : o.voicing;
  for (t = 0; t < T; t++){
    var pb = path[t], f0c = eng.fmin * Math.pow(2, pb * eng.binCents / 1200);
    var y0 = pb > 0 ? obs[t * nBins + pb - 1] : 0, y1 = obs[t * nBins + pb], y2 = pb < nBins - 1 ? obs[t * nBins + pb + 1] : 0;
    var dn2 = y0 - 2 * y1 + y2, dsh = dn2 !== 0 ? 0.5 * (y0 - y2) / dn2 : 0;
    if (dsh > 1 || dsh < -1) dsh = 0;
    f0c *= Math.pow(2, dsh * eng.binCents / 1200);
    var np = PN[t], pfS = PF.subarray(t * MP, t * MP + MP), paS = PA.subarray(t * MP, t * MP + MP);
    var f0r = np ? eng.refine(f0c, pfS, paS, np) : f0c;
    var salPath = (y1 / 255) * salMax[t];
    fracArr[t] = mass[t] > 0 ? salPath / mass[t] : 0; sp[t] = salPath;
    pitchAll[t] = hz2midi(f0r);
  }
  // hysteresis voicing: strong frames start a note, weaker frames only continue it while the pitch is contiguous
  var hiT = fracThr + 0.03, loT = Math.max(0.05, fracThr - 0.06), acc = new Uint8Array(T);
  var loud = function(i){ return maxA[i] >= 0.04 * p90; };
  for (t = 0; t < T; t++){
    if (loud(t) && fracArr[t] >= hiT) acc[t] = 1;
    else if (t > 0 && acc[t - 1] && loud(t) && fracArr[t] >= loT && Math.abs(pitchAll[t] - pitchAll[t - 1]) < 1.0) acc[t] = 1;
  }
  for (t = T - 2; t >= 0; t--){
    if (!acc[t] && acc[t + 1] && loud(t) && fracArr[t] >= loT && Math.abs(pitchAll[t] - pitchAll[t + 1]) < 1.0) acc[t] = 1;
  }
  for (t = 0; t < T; t++) midi[t] = acc[t] ? pitchAll[t] : NaN;
  yield 0.95;
  return { midi: midi, hop: hop / sr, t0: (N / 2) / sr, frac: fracArr, sp: sp };
}

/* ---------- SOLO (single melodic line) ---------- */
function* soloGen(x, sr, o){
  o = o || {};
  var W = 1024, hop = Math.round(0.010 * sr);
  var mpm = new MPM(sr, W, o.loHz || 70, o.hiHz || 1000, 0.9);
  var T = Math.max(1, Math.floor((x.length - W) / hop) + 1);
  var midi = new Float32Array(T), clar = new Float32Array(T), rmsA = new Float32Array(T);
  var last = now(), t;
  for (t = 0; t < T; t++){
    var r = mpm.detect(x, t * hop);
    rmsA[t] = r.rms; clar[t] = r.clarity;
    midi[t] = r.f0 > 0 ? hz2midi(r.f0) : NaN;
    if ((t & 63) === 0 && now() - last > 40){ last = now(); yield 0.95 * t / T; }
  }
  var sr2 = Array.prototype.slice.call(rmsA).sort(function(a, b){ return a - b; });
  var p90 = sr2[Math.floor(0.9 * (T - 1))] || 1e-9;
  var cthr = o.clarity === undefined ? 0.85 : o.clarity;
  var floor = Math.max(0.002, 0.04 * p90);
  for (t = 0; t < T; t++) if (!(clar[t] >= cthr && rmsA[t] >= floor)) midi[t] = NaN;
  yield 0.97;
  return { midi: midi, hop: hop / sr, t0: (W / 2) / sr, frac: clar };
}

/* ---------- track clean-up ---------- */
function cleanTrack(midi, hopSec, o){
  o = o || {};
  var T = midi.length, i, j;
  var out = new Float32Array(midi);
  // 1. split into continuous segments
  var segs = [], st = -1;
  for (i = 0; i < T; i++){
    var v = out[i];
    if (isNaN(v)){ if (st >= 0){ segs.push([st, i - 1]); st = -1; } continue; }
    if (st < 0) st = i;
    else if (Math.abs(v - out[i - 1]) > 3){ segs.push([st, i - 1]); st = i; }
  }
  if (st >= 0) segs.push([st, T - 1]);
  var med = segs.map(function(s){ return median(out.subarray(s[0], s[1] + 1)); });
  // 2. short segments sitting one octave (or two) away from a longer neighbour -> octave errors
  var maxGap = Math.round(0.3 / hopSec), shortLen = Math.round(0.08 / hopSec);
  for (i = 0; i < segs.length; i++){
    var len = segs[i][1] - segs[i][0] + 1;
    if (len >= shortLen) continue;
    var ref = NaN;
    for (j = i - 1; j >= 0 && j >= i - 3; j--){
      if (segs[i][0] - segs[j][1] > maxGap) break;
      if (segs[j][1] - segs[j][0] + 1 >= shortLen){ ref = med[j]; break; }
    }
    if (isNaN(ref)) for (j = i + 1; j < segs.length && j <= i + 3; j++){
      if (segs[j][0] - segs[i][1] > maxGap) break;
      if (segs[j][1] - segs[j][0] + 1 >= shortLen){ ref = med[j]; break; }
    }
    if (isNaN(ref)) continue;
    var d = med[i] - ref, k = Math.round(d / 12);
    if (k !== 0 && Math.abs(d - 12 * k) < 0.8){
      for (var q = segs[i][0]; q <= segs[i][1]; q++) out[q] -= 12 * k;
    }
  }
  // 3. 5-frame median on voiced runs
  var sm = new Float32Array(out), buf = [];
  for (i = 0; i < T; i++){
    if (isNaN(out[i])) continue;
    buf.length = 0;
    for (j = -2; j <= 2; j++){ var q2 = i + j; if (q2 >= 0 && q2 < T && !isNaN(out[q2]) && Math.abs(out[q2] - out[i]) < 3) buf.push(out[q2]); }
    sm[i] = median(buf);
  }
  out = sm;
  // 4. bridge very short gaps, drop very short voiced blips
  var gapMax = Math.round((o.bridge || 0.06) / hopSec), minRun = Math.round((o.minRun || 0.05) / hopSec);
  i = 0;
  while (i < T){
    if (isNaN(out[i])){
      var e = i; while (e < T && isNaN(out[e])) e++;
      if (i > 0 && e < T && e - i <= gapMax && Math.abs(out[e] - out[i - 1]) < 2){
        for (j = i; j < e; j++) out[j] = out[i - 1] + (out[e] - out[i - 1]) * (j - i + 1) / (e - i + 1);
 }
      i = e;
    } else i++;
  }
  i = 0;
  while (i < T){
    if (!isNaN(out[i])){
      var e2 = i; while (e2 < T && !isNaN(out[e2])) e2++;
      if (e2 - i < minRun) for (j = i; j < e2; j++) out[j] = NaN;
      i = e2;
    } else i++;
  }
  return out;
}

/* ---------- notes (swaras) from a pitch track and a Sa ---------- */
function buildNotes(midi, hopSec, t0, saMidi, o){
  o = o || {};
  var minNote = o.minNote === undefined ? 0.08 : o.minNote, hys = 0.15, T = midi.length;
  var ev = [], cur = null, i;
  function close(endIdx){
    if (!cur) return;
    cur.t1 = t0 + (endIdx + 1) * hopSec;
    cur.cents = cur.frames ? cur.sumC / cur.frames : 0;
    cur.dur = cur.t1 - cur.t0;
    ev.push(cur); cur = null;
  }
  for (i = 0; i < T; i++){
    var m = midi[i];
    if (isNaN(m)){ close(i - 1); continue; }
    var s = m - saMidi;
    if (!cur || Math.abs(s - cur.n) > 0.5 + hys){
      close(i - 1);
      var n = Math.round(s);
      cur = { n: n, idx: ((n % 12) + 12) % 12, oct: Math.floor(n / 12), t0: t0 + i * hopSec, t1: 0, frames: 0, sumC: 0, cents: 0, dur: 0, passing: false };
    }
    cur.frames++; cur.sumC += (s - cur.n) * 100;
  }
  close(T - 1);
  ev.forEach(function(e){ e.passing = e.dur < minNote; });
  // merge consecutive same-swara events once passing ornaments are set aside
  var main = [];
  ev.forEach(function(e){
    if (e.passing) return;
    var p = main[main.length - 1];
    if (p && p.n === e.n && e.t0 - p.t1 < 0.15){
      var tf = p.frames + e.frames;
      p.cents = (p.cents * p.frames + e.cents * e.frames) / tf; p.frames = tf; p.t1 = e.t1; p.dur = p.t1 - p.t0;
    } else main.push({ n: e.n, idx: e.idx, oct: e.oct, t0: e.t0, t1: e.t1, frames: e.frames, cents: e.cents, dur: e.dur, passing: false });
  });
  return { all: ev, main: main };
}

/* ---------- pitch-class histogram + Sa suggestion ---------- */
var SA_TEMPLATE = [1.0, 0.10, 0.55, 0.25, 0.50, 0.45, 0.15, 0.90, 0.15, 0.50, 0.25, 0.40];
function estimateSa(midi){
  var sx = 0, sy = 0, n = 0, i;
  for (i = 0; i < midi.length; i++){
    var m = midi[i]; if (isNaN(m)) continue;
    var ph = PI2 * (m - Math.floor(m)); sx += Math.cos(ph); sy += Math.sin(ph); n++;
  }
  if (n < 30) return null;
  var delta = Math.atan2(sy, sx) / PI2; // fractional-semitone offset of the singer's grid
  var P = new Float64Array(12), all = [];
  for (i = 0; i < midi.length; i++){
    var m2 = midi[i]; if (isNaN(m2)) continue;
    var pc = ((Math.round(m2 - delta) % 12) + 12) % 12; P[pc] += 1; all.push(m2);
  }
  var tot = 0; for (i = 0; i < 12; i++) tot += P[i];
  for (i = 0; i < 12; i++) P[i] /= tot;
  var scores = [];
  for (var c = 0; c < 12; c++){
    var sc = 0; for (var k = 0; k < 12; k++) sc += P[(c + k) % 12] * SA_TEMPLATE[k];
    scores.push(sc);
  }
  var order = scores.map(function(v, idx){ return idx; }).sort(function(a, b){ return scores[b] - scores[a]; });
  var best = order[0], conf = scores[order[0]] / (scores[order[1]] || 1e-9);
  var mid = median(all), base = best + delta;
  var saM = base + 12 * Math.floor((mid - base) / 12);
  if (mid - saM > 9) saM += 0; // keep highest Sa at or below the median pitch
  return { saMidi: saM, pc: best, delta: delta, confidence: conf, profile: P };
}

/* ---------- distribution of swaras (seconds per swara) ---------- */
function swaraHistogram(midi, hopSec, saMidi){
  var h = new Float64Array(12), tot = 0;
  for (var i = 0; i < midi.length; i++){
    if (isNaN(midi[i])) continue;
    var n = Math.round(midi[i] - saMidi), idx = ((n % 12) + 12) % 12;
    h[idx] += hopSec; tot += hopSec;
  }
  return { sec: h, total: tot };
}

return {
  FFT: FFT, MPM: MPM, SalienceEngine: SalienceEngine, decimate: decimate, hann: hann,
  melodyGen: melodyGen, soloGen: soloGen, cleanTrack: cleanTrack, buildNotes: buildNotes,
  estimateSa: estimateSa, swaraHistogram: swaraHistogram, rangePrior: rangePrior,
  hz2midi: hz2midi, midi2hz: midi2hz, median: median, nextPow2: nextPow2
};
})();
/* ==================== /DSP CORE ==================== */
    
