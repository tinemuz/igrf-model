/*
 * MIT License
 *
 * Copyright (c) 2025 tinemuz
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * This file contains derivative works. See the NOTICE file for attribution and
 * original license information:
 * https://github.com/tinemuz/igrf-model/blob/main/NOTICE
 */
package com.github.tinemuz;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * IGRF magnetic model evaluator.
 *
 * <p>This class evaluates the Earth's magnetic field using the published IGRF
 * spherical harmonic coefficients (degrees up to 13). It provides a single
 * main entry point, {@link #compute}, which returns a small immutable
 * {@link Field} object containing the field vector and common derived
 * quantities (declination, inclination, intensities).</p>
 *
 * <p>Inputs are standard WGS-84 geodetic latitude/longitude (degrees),
 * altitude in meters above mean sea level, and time as UTC epoch
 * milliseconds. Coefficients are loaded from the classpath resource
 * <code>igrfcoeffs.txt</code>. Call {@link #preload()} once at startup.</p>
 */
public final class IGRFModel {
    private static final Logger log = LoggerFactory.getLogger(IGRFModel.class);
    // WGS-84 ellipsoid parameters (km) and IAU reference radius (km)
    private static final float A_KM = 6378.137f;
    private static final float B_KM = 6356.7523142f;
    private static final float RE_KM = 6371.2f;
    private static final int MAX_N = 13;
    private static final float[][] SCHMIDT = computeSchmidtQuasiNormFactors(MAX_N + 1);
    private static final double MS_PER_YEAR_365 = 365.0 * 24.0 * 60.0 * 60.0 * 1000.0;
    
    /**
     * Per-thread scratch space so computations avoid frequent allocations.
     * Keeping a small workspace object on each thread improves speed for
     * repeated evaluations.
     */
    @SuppressWarnings("squid:S5164")
    private static final ThreadLocal<Workspace> WORK = ThreadLocal.withInitial(Workspace::new);
    private static volatile boolean loaded = false;
    private static double[] epochs; // available coefficient epochs (decimal years)
    private static float[][][] gCoeffs; // g[n][m][epochIndex]
    private static float[][][] hCoeffs; // h[n][m][epochIndex]
    private static float[][] svG; // secular variation for g (degrees <= 8)
    private static float[][] svH; // secular variation for h (degrees <= 8)
    private static float baseYearAnchor; // anchor year for decimal-year math
    private static long baseYearMillis; // millis for Jan 1 of anchor year
    private static volatile boolean warnedFutureBeyond5Years = false;

    private IGRFModel() {}

    /**
     * Evaluate the magnetic field at a WGS-84 geodetic location and time.
     *
     * Returns a {@link Field} with components in nanoteslas and angles in
     * degrees. The returned components are in the NED (north, east, down)
     * convention.
     *
     * This is the single public API most callers need.
     *
     * @param gdLatitudeDeg  geodetic latitude (degrees, north positive)
     * @param gdLongitudeDeg geodetic longitude (degrees, east positive)
     * @param altitudeMeters altitude above mean sea level (meters)
     * @param epochMillis    UTC time as epoch milliseconds
     * @return Field carrying X (north), Y (east), Z (down), H (horizontal),
     *         F (total), declination and inclination
     * @throws IllegalStateException if coefficient data cannot be loaded
     */
    public static Field compute(
            double gdLatitudeDeg, double gdLongitudeDeg, double altitudeMeters, long epochMillis) {
        ensureLoaded();
        
        // Keep latitude slightly inside [-90, 90] to avoid numeric issues at poles
        float latDeg = (float) Math.max(-90.0 + 1e-5, Math.min(gdLatitudeDeg, 90.0 - 1e-5));
        float lonDeg = (float) gdLongitudeDeg;
        
        // Convert geodetic coordinates (lat/lon/alt) to geocentric (radians, km)
        Geocentric gc = geodeticToGeocentric(latDeg, lonDeg, (float) (altitudeMeters / 1000.0));
        
        // Convert time to a decimal year anchored at the last epoch in the file
        float decimalYear =
                baseYearAnchor + (float) ((epochMillis - baseYearMillis) / MS_PER_YEAR_365);
        
        // Thread-local workspace (pre-allocated arrays we reuse)
        Workspace ws = WORK.get();

        // STEP 1: If colatitude changed, recompute Legendre functions and derivatives
        float theta = (float) (Math.PI / 2.0 - gc.latRad);
        if (theta != ws.lastTheta) {
            fillLegendre(ws, theta);
            ws.lastTheta = theta;
        }
        
        // STEP 2: If radius changed, recompute powers (a/r)^(n+1)
        if (gc.radiusKm != ws.lastRadiusKm) {
            ws.relPow[0] = 1.0f;
            ws.relPow[1] = RE_KM / gc.radiusKm;
            for (int i = 2; i < ws.relPow.length; i++)
                ws.relPow[i] = ws.relPow[i - 1] * ws.relPow[1];
            ws.lastRadiusKm = gc.radiusKm;
        }
        
        // STEP 3: If longitude changed, update sin(m*lon)/cos(m*lon) series
        if (gc.lonRad != ws.lastLonRad) {
            float s = (float) Math.sin(gc.lonRad);
            float c = (float) Math.cos(gc.lonRad);
            ws.sinMLon[0] = 0.0f;
            ws.cosMLon[0] = 1.0f;
            ws.sinMLon[1] = s;
            ws.cosMLon[1] = c;
            // Build the series using angle addition to avoid repeated trig calls
            for (int m = 2; m <= MAX_N; m++) {
                float sp = ws.sinMLon[m - 1];
                float cp = ws.cosMLon[m - 1];
                ws.sinMLon[m] = sp * c + cp * s;
                ws.cosMLon[m] = cp * c - sp * s;
            }
            ws.lastLonRad = gc.lonRad;
        }
        
        // Precompute 1/cos(lat) used for east component (guarded by lat clamp above)
        float invCosLat = 1.0f / (float) Math.cos(gc.latRad);
        final int E = epochs.length;
        final double lastEpoch = epochs[E - 1];
        
        // Accumulate geocentric components (these are gradients of the potential)
        float gcX = 0f;  // component along colatitude (θ)
        float gcY = 0f;  // component along longitude (φ)
        float gcZ = 0f;  // radial component (r)

        // STEP 4: Evaluate spherical harmonic sum. Coefficients depend on time.

        // Case A: Requested time is before the earliest epoch — use first column
        if (decimalYear <= epochs[0]) {
            for (int n = 1; n <= MAX_N; n++) {
                float rn = ws.relPow[n + 2];  // (a/r)^(n+1)
                float[] pS = ws.ps[n];        // Schmidt-normalized P_n^m
                float[] dpS = ws.pds[n];      // derivative dP_n^m/dθ
                for (int m = 0; m <= n; m++) {
                    float gnm = gCoeffs[n][m][0];
                    float hnm = hCoeffs[n][m][0];
                    float a = gnm * ws.cosMLon[m] + hnm * ws.sinMLon[m];
                    float b = gnm * ws.sinMLon[m] - hnm * ws.cosMLon[m];
                    gcX += rn * a * dpS[m];
                    if (m != 0) gcY += rn * m * b * pS[m] * invCosLat;
                    gcZ -= rn * (n + 1) * a * pS[m];
                }
            }
        } else if (decimalYear <= lastEpoch) {
            // Case B: Between epochs — linearly interpolate coefficients between nearest epochs
            int hi = upperBound(epochs, decimalYear);
            int lo = hi - 1;
            float t = (float) ((decimalYear - epochs[lo]) / (epochs[hi] - epochs[lo]));
            for (int n = 1; n <= MAX_N; n++) {
                float rn = ws.relPow[n + 2];
                float[] pS = ws.ps[n];
                float[] dpS = ws.pds[n];
                for (int m = 0; m <= n; m++) {
                    float g0 = gCoeffs[n][m][lo];
                    float g1 = gCoeffs[n][m][hi];
                    float h0 = hCoeffs[n][m][lo];
                    float h1 = hCoeffs[n][m][hi];
                    float gnm = g0 + t * (g1 - g0);
                    float hnm = h0 + t * (h1 - h0);
                    float a = gnm * ws.cosMLon[m] + hnm * ws.sinMLon[m];
                    float b = gnm * ws.sinMLon[m] - hnm * ws.cosMLon[m];
                    gcX += rn * a * dpS[m];
                    if (m != 0) gcY += rn * m * b * pS[m] * invCosLat;
                    gcZ -= rn * (n + 1) * a * pS[m];
                }
            }
        } else {
            // Case C: After the last epoch — use last epoch plus secular variation (SV)
            // SV is only published for degrees up to 8, and we clamp extrapolation to 5 years
            float yearsAhead = (float) (decimalYear - lastEpoch);
            float dt = yearsAhead > 5.0f ? 5.0f : yearsAhead;  // clamp to max +5 years
            if (yearsAhead > 5.0f && !warnedFutureBeyond5Years) {
                synchronized (IGRFModel.class) {
                    if (!warnedFutureBeyond5Years) {
                        warnedFutureBeyond5Years = true;
                        log.warn(
                                "Requested time is {} years beyond latest epoch {}; "
                                    + "results are clamped to {}. Consider updating igrfcoeffs.txt",
                                String.format("%.2f", yearsAhead), 
                                String.format("%.1f", lastEpoch), 
                                String.format("%.1f", lastEpoch + 5.0));
                    }
                }
            }
            for (int n = 1; n <= MAX_N; n++) {
                float rn = ws.relPow[n + 2];
                float[] pS = ws.ps[n];
                float[] dpS = ws.pds[n];
                boolean useSv = n <= 8; // SV only available for lower degrees
                for (int m = 0; m <= n; m++) {
                    float baseG = gCoeffs[n][m][E - 1];
                    float baseH = hCoeffs[n][m][E - 1];
                    float gnm = useSv ? baseG + dt * svG[n][m] : baseG;
                    float hnm = useSv ? baseH + dt * svH[n][m] : baseH;
                    float a = gnm * ws.cosMLon[m] + hnm * ws.sinMLon[m];
                    float b = gnm * ws.sinMLon[m] - hnm * ws.cosMLon[m];
                    gcX += rn * a * dpS[m];
                    if (m != 0) gcY += rn * m * b * pS[m] * invCosLat;
                    gcZ -= rn * (n + 1) * a * pS[m];
                }
            }
        }
        
        // STEP 5: Rotate geocentric components to geodetic NED frame.
        // The geodetic latitude differs slightly from geocentric latitude;
        // apply that small rotation before returning north/east/down.
        float latDiff = (float) (Math.toRadians(latDeg) - gc.latRad);
        float cosLd = (float) Math.cos(latDiff);
        float sinLd = (float) Math.sin(latDiff);
        // X_ned = X_gc*cos(Δlat) + Z_gc*sin(Δlat)
        // Y_ned = Y_gc
        // Z_ned = -X_gc*sin(Δlat) + Z_gc*cos(Δlat)
        return new Field(gcX * cosLd + gcZ * sinLd, gcY, -gcX * sinLd + gcZ * cosLd);
    }

    /**
     * Ensure coefficient arrays are loaded. Safe to call repeatedly; the
     *      * first caller will load data from the classpath. This method is synchronized
     *      * to avoid race conditions during initialization.
     */
    private static synchronized void ensureLoaded() {
        if (loaded) return;
        loadCoefficientsFromResource();
        loaded = true;
    }

    /**
     * Load coefficient data and prepare internal arrays. Call this at startup
     * if you want to detect missing or invalid coefficient files early.
     */
    public static void preload() {
        ensureLoaded();
        if (!warnedFutureBeyond5Years && epochs != null && epochs.length > 0) {
            double lastEpoch = epochs[epochs.length - 1];
            long now = System.currentTimeMillis();
            float nowYear = baseYearAnchor + (float) ((now - baseYearMillis) / MS_PER_YEAR_365);
            float yearsAhead = (float) (nowYear - lastEpoch);
            if (yearsAhead > 5.0f) {
                synchronized (IGRFModel.class) {
                    if (!warnedFutureBeyond5Years) {
                        warnedFutureBeyond5Years = true;
                        log.warn(
                                "Current time is {} years beyond latest epoch {}; "
                                    + "results will clamp to {}. Consider updating igrfcoeffs.txt",
                                String.format("%.2f", yearsAhead), 
                                String.format("%.1f", lastEpoch), 
                                String.format("%.1f", lastEpoch + 5.0));
                    }
                }
            }
        }
    }

    /**
     * Parse a token from the header line and add it to the years list if it's
     * a numeric epoch. Known non-year tokens like "SV" or "g/h" are ignored.
     */
    private static void parseYearToken(String token, List<Double> years) {
        if (token.equalsIgnoreCase("SV") || token.equalsIgnoreCase("g/h")) {
            return;
        }
        try {
            years.add(Double.parseDouble(token));
        } catch (NumberFormatException e) {
            // Ignore tokens that aren't parseable as a year — header robustness
        }
    }

    /**
     * Load IGRF coefficient table from the classpath resource `igrfcoeffs.txt`.
     * The file format is the standard table of g/h rows with columns for each
     * epoch; the method fills gCoeffs, hCoeffs and optional secular variation
     * arrays. Any problems reading or parsing the file will throw an
     * IllegalStateException.
     */
    private static void loadCoefficientsFromResource() {
        InputStream in = IGRFModel.class.getClassLoader().getResourceAsStream("igrfcoeffs.txt");
        if (in == null) {
            log.error("IGRF coefficients file 'igrfcoeffs.txt' not found on classpath");
            throw new IllegalStateException(
                    "IGRF coefficients file 'igrfcoeffs.txt' not found on classpath");
        }
        try (BufferedReader br = new BufferedReader(new InputStreamReader(in))) {
            List<Double> years = new ArrayList<>();
            List<float[]> gRows = new ArrayList<>();
            List<float[]> hRows = new ArrayList<>();
            boolean headerParsed = false;
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) continue; // skip blanks/comments
                if (!headerParsed && line.startsWith("g/h")) {
                    String[] toks = line.split("\\s+");
                    years.clear();
                    for (String t : toks) {
                        parseYearToken(t, years);
                    }
                    headerParsed = !years.isEmpty();
                    continue;
                }
                char kind = line.charAt(0);
                if (kind != 'g' && kind != 'h') continue; // only care about g/h data rows
                String[] toks = line.split("\\s+");
                if (toks.length < 4) continue; // malformed row
                int n = Integer.parseInt(toks[1]);
                int m = Integer.parseInt(toks[2]);
                int startIdx = 3;
                int count = toks.length - startIdx;
                float[] row = new float[3 + count];
                row[0] = n;
                row[1] = m;
                row[2] = (kind == 'g') ? 1f : 2f; // marker: 1 for g, 2 for h
                for (int i = 0; i < count; i++) row[3 + i] = Float.parseFloat(toks[startIdx + i]);
                if (kind == 'g') gRows.add(row);
                else hRows.add(row);
            }
            if (years.isEmpty()) {
                // If header had no years, infer them from row length using a simple 5-year step
                int dataEpochs = gRows.isEmpty() ? 0 : gRows.get(0).length - 3;
                List<Double> defaultYears = new ArrayList<>();
                for (int i = 0; i < dataEpochs; i++) {
                    defaultYears.add(1900.0 + i * 5.0);
                }
                years = defaultYears;
            }
            epochs = years.stream().mapToDouble(Double::doubleValue).toArray();
            int epochCount = epochs.length;
            int latestYear = (int) Math.round(epochs[epochCount - 1]);
            baseYearAnchor = (float) epochs[epochCount - 1];
            baseYearMillis = millisAtUtcJan1(latestYear);
            gCoeffs = new float[MAX_N + 1][MAX_N + 1][epochCount];
            hCoeffs = new float[MAX_N + 1][MAX_N + 1][epochCount];
            svG = new float[MAX_N + 1][MAX_N + 1];
            svH = new float[MAX_N + 1][MAX_N + 1];
            for (float[] row : gRows) {
                int n = (int) row[0];
                int m = (int) row[1];
                int evals = row.length - 3;
                System.arraycopy(row, 3, gCoeffs[n][m], 0, epochCount);
                if (evals > epochCount) svG[n][m] = row[3 + epochCount];
            }
            for (float[] row : hRows) {
                int n = (int) row[0];
                int m = (int) row[1];
                int evals = row.length - 3;
                System.arraycopy(row, 3, hCoeffs[n][m], 0, epochCount);
                if (evals > epochCount) svH[n][m] = row[3 + epochCount];
            }
        } catch (IOException e) {
            log.error("Failed to read IGRF coefficients file", e);
            throw new IllegalStateException("Failed to read IGRF coefficients file", e);
        } catch (Exception e) {
            log.error("Failed to parse IGRF coefficients file", e);
            throw new IllegalStateException("Failed to parse IGRF coefficients file", e);
        }
    }

    /**
     * Milliseconds at UTC 00:00:00 on January 1 of the given year.
     * Uses a proleptic Gregorian calendar calculation.
     */
    private static long millisAtUtcJan1(int year) {
        return daysFrom1970ToYearStart(year) * 86_400_000L;
    }

    /**
     * Days from 1970-01-01 to the start of the given year.
     *
     * This uses a simple arithmetic calculation (365-day years plus leap days)
     * and is intended for coarse epoch math (millisecond resolution for
     * year boundaries) used by the model anchoring code.
     */
    private static long daysFrom1970ToYearStart(int year) {
        if (year == 1970) return 0L;
        long diff = year - 1970L;
        return 365L * diff + (leapsBefore(year) - leapsBefore(1970));
    }

    /**
     * Count leap days strictly before the start of the given year.
     *
     * Uses the proleptic Gregorian rules: a year is a leap year if
     * divisible by 4, except centuries not divisible by 400.
     */
    private static long leapsBefore(int year) {
        long y = year - 1L;
        return y / 4L - y / 100L + y / 400L;
    }

    /**
     * Find the first index in a sorted array where arr[index] >= x.
     * If x is larger than all entries, returns the last index.
     */
    private static int upperBound(double[] arr, double x) {
        int lo = 0;
        int hi = arr.length - 1;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (arr[mid] < x) lo = mid + 1;
            else hi = mid;
        }
        return lo;
    }

    /**
     * Convert geodetic (WGS-84) lat/lon/alt into geocentric latitude (radians),
     * longitude (radians), and radius (km). Altitude is expected in km.
     */
    private static Geocentric geodeticToGeocentric(float gdLatDeg, float gdLonDeg, float altKm) {
        float a2 = A_KM * A_KM;
        float b2 = B_KM * B_KM;
        double gdLatRad = Math.toRadians(gdLatDeg);
        float clat = (float) Math.cos(gdLatRad);
        float slat = (float) Math.sin(gdLatRad);
        float tlat = slat / clat;
        float rho = (float) Math.sqrt(a2 * clat * clat + b2 * slat * slat);
        float latRad = (float) Math.atan(tlat * (rho * altKm + b2) / (rho * altKm + a2));
        float lonRad = (float) Math.toRadians(gdLonDeg);
        float radSq =
                altKm * altKm
                        + 2 * altKm * rho
                        + (a2 * a2 * clat * clat + b2 * b2 * slat * slat)
                                / (a2 * clat * clat + b2 * slat * slat);
        float radiusKm = (float) Math.sqrt(radSq);
        return new Geocentric(latRad, lonRad, radiusKm);
    }

    /**
     * Compute associated Legendre functions P_n^m(cos(theta)) and their
     * theta-derivatives, then apply Schmidt quasi-normalization. Results are
     * placed in the {@link Workspace} arrays used by the main evaluate loop.
     */
    private static void fillLegendre(Workspace ws, float thetaRad) {
        float cos = (float) Math.cos(thetaRad);
        float sin = (float) Math.sin(thetaRad);
        ws.p[0][0] = 1.0f;
        ws.pDeriv[0][0] = 0.0f;
        for (int n = 1; n <= MAX_N; n++) {
            for (int m = 0; m <= n; m++) {
                if (n == m) {
                    ws.p[n][m] = sin * ws.p[n - 1][m - 1];
                    ws.pDeriv[n][m] = cos * ws.p[n - 1][m - 1] + sin * ws.pDeriv[n - 1][m - 1];
                } else if (n == 1 || m == n - 1) {
                    ws.p[n][m] = cos * ws.p[n - 1][m];
                    ws.pDeriv[n][m] = -sin * ws.p[n - 1][m] + cos * ws.pDeriv[n - 1][m];
                } else {
                    float k = ((n - 1f) * (n - 1f) - m * m) / ((2f * n - 1f) * (2f * n - 3f));
                    ws.p[n][m] = cos * ws.p[n - 1][m] - k * ws.p[n - 2][m];
                    ws.pDeriv[n][m] =
                            -sin * ws.p[n - 1][m]
                                    + cos * ws.pDeriv[n - 1][m]
                                    - k * ws.pDeriv[n - 2][m];
                }
            }
        }
        for (int n = 0; n <= MAX_N; n++) {
            float[] sch = SCHMIDT[n];
            for (int m = 0; m <= n; m++) {
                ws.ps[n][m] = ws.p[n][m] * sch[m];
                ws.pds[n][m] = ws.pDeriv[n][m] * sch[m];
            }
        }
    }

    /**
     * Compute Schmidt quasi-normalization factors used to scale the associated
     * Legendre functions. The input is inclusive maximum order (e.g. MAX_N+1).
     */
    private static float[][] computeSchmidtQuasiNormFactors(int maxNInclusive) {
        float[][] s = new float[maxNInclusive][];
        s[0] = new float[] {1.0f};
        for (int n = 1; n < maxNInclusive; n++) {
            s[n] = new float[n + 1];
            s[n][0] = s[n - 1][0] * (2f * n - 1f) / n;
            for (int m = 1; m <= n; m++)
                s[n][m] =
                        s[n][m - 1]
                                * (float) Math.sqrt((n - m + 1f) * (m == 1 ? 2f : 1f) / (n + m));
        }
        return s;
    }

    /**
     * Small immutable result object returned by {@link #compute}.
     * All linear quantities are in nanoteslas (nT); angles are degrees.
     */
    public static final class Field {
        /** X component (north) in nanoteslas (nT). */
        public final float xNorthNt;

        /** Y component (east) in nanoteslas (nT). */
        public final float yEastNt;

        /** Z component (down) in nanoteslas (nT). Positive downward. */
        public final float zDownNt;

        /** Horizontal intensity H = sqrt(X^2 + Y^2) in nanoteslas (nT). */
        public final float hHorizontalNt;

        /** Total intensity F = sqrt(X^2 + Y^2 + Z^2) in nanoteslas (nT). */
        public final float fTotalNt;

        /** Declination D = atan2(Y, X) in degrees, east of geographic north. */
        public final float declinationDeg;

        /** Inclination I (dip) = atan2(Z, H) in degrees, positive downward. */
        public final float inclinationDeg;

        private Field(float x, float y, float z) {
            this.xNorthNt = x;
            this.yEastNt = y;
            this.zDownNt = z;
            this.hHorizontalNt = (float) Math.hypot(x, y);
            this.fTotalNt = (float) Math.sqrt(x * x + y * y + z * z);
            this.declinationDeg = (float) Math.toDegrees(Math.atan2(y, x));
            this.inclinationDeg = (float) Math.toDegrees(Math.atan2(z, this.hHorizontalNt));
        }
    }

    // Simple holder for geocentric lat/lon/radius
    private record Geocentric(float latRad, float lonRad, float radiusKm) {}

    // Workspace: pre-allocated arrays used during a single evaluation
    private static final class Workspace {
        final float[][] p = new float[MAX_N + 1][MAX_N + 1];
        final float[][] pDeriv = new float[MAX_N + 1][MAX_N + 1];
        final float[][] ps = new float[MAX_N + 1][MAX_N + 1];
        final float[][] pds = new float[MAX_N + 1][MAX_N + 1];
        final float[] sinMLon = new float[MAX_N + 1];
        final float[] cosMLon = new float[MAX_N + 1];
        final float[] relPow = new float[MAX_N + 3];
        float lastTheta = Float.NaN;
        float lastLonRad = Float.NaN;
        float lastRadiusKm = Float.NaN;
    }
}
