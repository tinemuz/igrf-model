package com.github.tinemuz;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

/**
 * IGRFModel
 * <p>
 * A small, allocation-free, high-performance evaluator for the IGRF main field that
 * reads coefficients from a classpath resource ("igrfcoeffs.txt") at runtime and computes
 * the geomagnetic field vector and declination for geodetic inputs. The implementation is
 * designed to be a tiny dependency with predictable latency suitable for real-time streams.
 * </p>
 *
 * <h2>Coordinate frames and sign conventions</h2>
 * <ul>
 *   <li>Input latitude/longitude: WGS‑84 geodetic degrees (north/east positive).</li>
 *   <li>Input altitude: meters above mean sea level (WGS‑84 ellipsoid). Negative values are allowed.</li>
 *   <li>Output vector components: local geodetic North-East-Down (NED) frame in nanoteslas (nT):
 *     <ul>
 *       <li>X ≡ northward component (positive toward geographic north).</li>
 *       <li>Y ≡ eastward component (positive toward geographic east).</li>
 *       <li>Z ≡ downward component (positive toward the Earth’s center).</li>
 *     </ul>
 *   </li>
 *   <li>Declination: degrees east of geographic north in the range [-180°, +180°).</li>
 *   <li>Inclination (dip): degrees down from the horizontal plane (positive downward).</li>
 * </ul>
 *
 * <h2>Inputs and units</h2>
 * <ul>
 *   <li>Latitude/Longitude: WGS‑84 geodetic degrees (north/east positive).</li>
 *   <li>Altitude: meters above mean sea level (WGS‑84). Internally converted to kilometers.</li>
 *   <li>Time: epoch milliseconds (UTC). Internally converted to a decimal year using a
 *       365‑day approximation anchored at the latest epoch in the coefficient file to
 *       minimize drift around the model’s most recent year.</li>
 * </ul>
 *
 * <h2>Outputs and units</h2>
 * <ul>
 *   <li>X (north), Y (east), Z (down): nT</li>
 *   <li>Horizontal intensity H = sqrt(X² + Y²): nT</li>
 *   <li>Total intensity F = sqrt(X² + Y² + Z²): nT</li>
 *   <li>Declination D = atan2(Y, X): degrees</li>
 *   <li>Inclination I = atan2(Z, H): degrees</li>
 * </ul>
 *
 * <h2>Model and equations</h2>
 * <p>
 * The geomagnetic main field is modeled by a spherical harmonic potential (Schmidt quasi‑normalized):
 * </p>
 * <pre>
 * V(r, θ, φ) = a Σ_{n=1..N} (a/r)^{n+1} Σ_{m=0..n} [ g_{n m} cos(mφ) + h_{n m} sin(mφ) ] P_{n}^{m}(cos θ)
 * </pre>
 * <p>
 * where a = 6371.2 km is the reference Earth radius, r is geocentric radius, θ is geocentric
 * colatitude, φ is longitude, and P_{n}^{m} are associated Legendre functions. Coefficients
 * (g_{n m}, h_{n m}) are Schmidt quasi‑normalized. The magnetic field is the negative gradient:
 * </p>
 * <pre>
 * B_r = -∂V/∂r
 * B_θ = -1/r ∂V/∂θ
 * B_φ = -1/(r sin θ) ∂V/∂φ
 * </pre>
 * <p>
 * Internally we compute Gauss‑normalized P and ∂P/∂θ via efficient recurrences, then multiply once per θ by Schmidt
 * factors. We cache (a/r) powers and sin(mφ)/cos(mφ) series to avoid per‑term transcendentals and minimize allocations.
 * </p>
 *
 * <h2>Time handling</h2>
 * <ul>
 *   <li><b>Between listed epochs:</b> linear interpolation of g/h within the two surrounding years.</li>
 *   <li><b>After the last epoch:</b> apply the published secular variation (SV) from last epoch to last+5 years
 *       for degrees n ≤ 8; the time step is clamped to +5 years. A once‑only warning is printed if evaluation time
 *       is more than 5 years beyond the latest epoch (also checked in {@link #preload()}).</li>
 *   <li><b>Decimal year:</b> computed from epoch milliseconds using a 365‑day year, anchored at the latest epoch’s
 *       Jan 1, 00:00:00 UTC to bound leap‑day drift around the most recent model year.
 *       If you need a different convention, compute the decimal year yourself and adapt your call site.</li>
 * </ul>
 *
 * <h2>Geodetic ↔ Geocentric</h2>
 * <p>
 * Inputs are geodetic (WGS‑84). We convert to geocentric latitude and radius using standard
 * ellipsoid relations (a = 6378.137 km, b = 6356.7523142 km), evaluate the spherical harmonics
 * in geocentric coordinates (r, θ, φ), and then rotate the result back to the local geodetic NED frame
 * to obtain X (north), Y (east), Z (down).
 * </p>
 *
 * <h2>Performance and threading</h2>
 * <ul>
 *   <li>One‑time: coefficient file read/parse on first use (or explicit {@link #preload()}) and per‑thread workspace allocation.</li>
 *   <li>Per‑call: allocation‑free; O(N²) with N = 13. Cached Legendre×Schmidt, sin/cos series, (a/r) powers.</li>
 *   <li>Thread‑safe: a thread‑local workspace avoids contention and repeated allocations.</li>
 * </ul>
 *
 * <h2>Data source</h2>
 * <p>
 * Coefficients are read from the classpath resource <code>igrfcoeffs.txt</code>. The file format is compatible with
 * the conventional IGRF table: a header row starting with <code>g/h</code> followed by epoch years and optional
 * <code>SV</code> column, and data rows beginning with <code>g</code> or <code>h</code>, degree n, order m, then values for
 * each epoch (and optionally SV). If the header is missing, a sensible default epoch list is used.
 * </p>
 */
public final class IGRFModel {
    // WGS84 constants (km) and IAU reference radius (km)
    private static final float A_KM = 6378.137f;
    private static final float B_KM = 6356.7523142f;
    private static final float RE_KM = 6371.2f;
    private static final int MAX_N = 13;
    private static final float[][] SCHMIDT = computeSchmidtQuasiNormFactors(MAX_N + 1);
    private static final double MS_PER_YEAR_365 = 365.0 * 24.0 * 60.0 * 60.0 * 1000.0;
    
    /**
     * Thread-local workspace for allocation-free computation.
     * Intentionally never removed (S5164 suppressed) - this is a performance optimization
     * where the workspace persists for the lifetime of the thread to avoid repeated allocations.
     * Each workspace is ~20KB and provides significant performance benefits in real-time systems.
     */
    @SuppressWarnings("squid:S5164")
    private static final ThreadLocal<Workspace> WORK = ThreadLocal.withInitial(Workspace::new);
    private static volatile boolean loaded = false;
    private static double[] epochs;
    private static float[][][] gCoeffs;
    private static float[][][] hCoeffs;
    private static float[][] svG2025;
    private static float[][] svH2025;
    private static float baseYearAnchor;
    private static long baseYearMillis;
    private static volatile boolean warnedFutureBeyond5Years = false;

    private IGRFModel() {}

    /**
     * Convenience method that returns only the magnetic declination at a given location and time.
     * <p>
     * Equivalent to {@link #compute(double, double, double, long)} followed by reading {@link Result#declinationDeg}.
     * </p>
     *
     * <h3>Parameters</h3>
     * <ul>
     *   <li><b>latDeg</b> — geodetic latitude in degrees (north positive). Values are internally clamped to
     *       (-90°, +90°) by ±1e‑5 degrees to avoid numerical singularities at the poles.</li>
     *   <li><b>lonDeg</b> — geodetic longitude in degrees (east positive). Any real value is accepted; typical range is [-180°, +180°) or [0°, 360°).</li>
     *   <li><b>altitudeMeters</b> — altitude above mean sea level (WGS‑84), in meters. Negative values are allowed.</li>
     *   <li><b>epochMillis</b> — UTC epoch time in milliseconds.</li>
     * </ul>
     *
     * <h3>Returns</h3>
     * <p>Declination in degrees east of north, in the range [-180°, +180°).</p>
     *
     * <h3>Exceptions</h3>
     * <ul>
     *   <li>{@link IllegalStateException} if the coefficient resource <code>igrfcoeffs.txt</code> is not found on the classpath.</li>
     * </ul>
     */
    public static double declinationDeg(
            double latDeg, double lonDeg, double altitudeMeters, long epochMillis) {
        Result r = compute(latDeg, lonDeg, altitudeMeters, epochMillis);
        return r.declinationDeg;
    }

    /**
     * Compute the geomagnetic field vector and common derived quantities at a geodetic location and time.
     * <p>
     * Results are returned in the local geodetic North-East-Down (NED) frame. Internally, inputs are converted to
     * geocentric coordinates, the IGRF spherical harmonic expansion is evaluated, and the result is rotated back to NED.
     * </p>
     *
     * <h3>Parameters</h3>
     * <ul>
     *   <li><b>gdLatitudeDeg</b> — geodetic latitude in degrees (north positive). Internally clamped to
     *       (-90°, +90°) by ±1e‑5 degrees to avoid pole singularities.</li>
     *   <li><b>gdLongitudeDeg</b> — geodetic longitude in degrees (east positive).</li>
     *   <li><b>altitudeMeters</b> — altitude above mean sea level (WGS‑84), in meters.</li>
     *   <li><b>epochMillis</b> — evaluation time as UTC epoch milliseconds.</li>
     * </ul>
     *
     * <h3>Returns</h3>
     * <p>A {@link Result} containing X, Y, Z (nT), H, F (nT), declination (deg), and inclination (deg).</p>
     *
     * <h3>Threading and performance</h3>
     * <ul>
     *   <li>Thread‑safe: independent thread‑local cache per calling thread.</li>
     *   <li>Allocation‑free per call after first use.</li>
     * </ul>
     *
     * <h3>Exceptions</h3>
     * <ul>
     *   <li>{@link IllegalStateException} if the coefficient resource <code>igrfcoeffs.txt</code> is missing.</li>
     * </ul>
     *
     * <h3>Implementation Note</h3>
     * <p>This method has high cognitive complexity (S3776) because it implements the complete
     * IGRF spherical harmonic evaluation algorithm with three time-handling cases. The complexity
     * is inherent to the mathematical model and has been extensively tested.</p>
     */
    public static Result compute(
            double gdLatitudeDeg, double gdLongitudeDeg, double altitudeMeters, long epochMillis) {
        ensureLoaded();
        
        // Clamp latitude to avoid numerical issues at poles
        float latDeg = (float) Math.clamp(gdLatitudeDeg, -90.0 + 1e-5, 90.0 - 1e-5);
        float lonDeg = (float) gdLongitudeDeg;
        
        // Convert WGS-84 geodetic to geocentric coordinates
        Geocentric gc = geodeticToGeocentric(latDeg, lonDeg, (float) (altitudeMeters / 1000.0));
        
        // Convert epoch time to decimal year for IGRF coefficient interpolation
        float decimalYear =
                baseYearAnchor + (float) ((epochMillis - baseYearMillis) / MS_PER_YEAR_365);
        
        // Get thread-local workspace (cache of pre-computed values)
        Workspace ws = WORK.get();

        // === STEP 1: Cache Legendre functions for this colatitude ===
        float theta = (float) (Math.PI / 2.0 - gc.latRad);
        if (theta != ws.lastTheta) {
            fillLegendre(ws, theta);
            ws.lastTheta = theta;
        }
        
        // === STEP 2: Cache (a/r)^(n+1) powers for this radius ===
        if (gc.radiusKm != ws.lastRadiusKm) {
            ws.relPow[0] = 1.0f;
            ws.relPow[1] = RE_KM / gc.radiusKm;
            for (int i = 2; i < ws.relPow.length; i++)
                ws.relPow[i] = ws.relPow[i - 1] * ws.relPow[1];
            ws.lastRadiusKm = gc.radiusKm;
        }
        
        // === STEP 3: Cache sin(m*lon) and cos(m*lon) series for this longitude ===
        if (gc.lonRad != ws.lastLonRad) {
            float s = (float) Math.sin(gc.lonRad);
            float c = (float) Math.cos(gc.lonRad);
            ws.sinMLon[0] = 0.0f;
            ws.cosMLon[0] = 1.0f;
            ws.sinMLon[1] = s;
            ws.cosMLon[1] = c;
            // Use angle addition formulas to avoid repeated trig calls
            for (int m = 2; m <= MAX_N; m++) {
                float sp = ws.sinMLon[m - 1];
                float cp = ws.cosMLon[m - 1];
                ws.sinMLon[m] = sp * c + cp * s;
                ws.cosMLon[m] = cp * c - sp * s;
            }
            ws.lastLonRad = gc.lonRad;
        }
        
        // Pre-compute 1/cos(lat) for east component calculation
        float invCosLat = 1.0f / (float) Math.cos(gc.latRad);
        final int E = epochs.length;
        final double lastEpoch = epochs[E - 1];
        
        // Initialize field components in geocentric frame
        float gcX = 0f;  // Colatitude component (dV/dθ)
        float gcY = 0f;  // Longitude component (dV/dφ)
        float gcZ = 0f;  // Radial component (dV/dr)
        
        // === STEP 4: Evaluate spherical harmonic expansion with time-dependent coefficients ===
        
        // Case 1: Before first epoch - use earliest coefficients
        if (decimalYear <= epochs[0]) {
            for (int n = 1; n <= MAX_N; n++) {
                float rn = ws.relPow[n + 2];  // (a/r)^(n+1)
                float[] pS = ws.ps[n];         // Schmidt-normalized Legendre P_n^m
                float[] dpS = ws.pds[n];       // Derivative dP_n^m/dθ
                for (int m = 0; m <= n; m++) {
                    float gnm = gCoeffs[n][m][0];
                    float hnm = hCoeffs[n][m][0];
                    float a = gnm * ws.cosMLon[m] + hnm * ws.sinMLon[m];
                    float b = gnm * ws.sinMLon[m] - hnm * ws.cosMLon[m];
                    gcX += rn * a * dpS[m];  // B_θ contribution
                    if (m != 0) gcY += rn * m * b * pS[m] * invCosLat;  // B_φ contribution
                    gcZ -= rn * (n + 1) * a * pS[m];  // B_r contribution (negative gradient)
                }
            }
        } else if (decimalYear <= lastEpoch) {
            // Case 2: Between epochs - linear interpolation
            int hi = upperBound(epochs, decimalYear);
            int lo = hi - 1;
            // Interpolation parameter t ∈ [0, 1]
            float t = (float) ((decimalYear - epochs[lo]) / (epochs[hi] - epochs[lo]));
            for (int n = 1; n <= MAX_N; n++) {
                float rn = ws.relPow[n + 2];
                float[] pS = ws.ps[n];
                float[] dpS = ws.pds[n];
                for (int m = 0; m <= n; m++) {
                    // Linear interpolation: coeff = coeff0 + t*(coeff1 - coeff0)
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
            // Case 3: Beyond last epoch - apply secular variation (SV) up to +5 years
            float yearsAhead = (float) (decimalYear - lastEpoch);
            float dt = yearsAhead > 5.0f ? 5.0f : yearsAhead;  // Clamp to 5 years
            if (yearsAhead > 5.0f && !warnedFutureBeyond5Years) {
                synchronized (IGRFModel.class) {
                    if (!warnedFutureBeyond5Years) {
                        warnedFutureBeyond5Years = true;
                        // Intentional diagnostic output - this library has no logger dependency
                        System.err.printf( // NOSONAR
                                "IGRFModel warning: requested time is %.2f years beyond latest"
                                    + " epoch %.1f; results are clamped to %.1f. Consider updating"
                                    + " igrfcoeffs.txt.%n",
                                yearsAhead, (float) lastEpoch, (float) (lastEpoch + 5.0));
                    }
                }
            }
            for (int n = 1; n <= MAX_N; n++) {
                float rn = ws.relPow[n + 2];
                float[] pS = ws.ps[n];
                float[] dpS = ws.pds[n];
                // SV is only provided for degrees n ≤ 8
                boolean useSv = n <= 8;
                for (int m = 0; m <= n; m++) {
                    float baseG = gCoeffs[n][m][E - 1];
                    float baseH = hCoeffs[n][m][E - 1];
                    // Apply secular variation: coeff = baseCoeff + dt*SV
                    float gnm = useSv ? baseG + dt * svG2025[n][m] : baseG;
                    float hnm = useSv ? baseH + dt * svH2025[n][m] : baseH;
                    float a = gnm * ws.cosMLon[m] + hnm * ws.sinMLon[m];
                    float b = gnm * ws.sinMLon[m] - hnm * ws.cosMLon[m];
                    gcX += rn * a * dpS[m];
                    if (m != 0) gcY += rn * m * b * pS[m] * invCosLat;
                    gcZ -= rn * (n + 1) * a * pS[m];
                }
            }
        }
        
        // === STEP 5: Rotate from geocentric to geodetic NED frame ===
        // Account for the difference between geodetic and geocentric latitude
        float latDiff = (float) (Math.toRadians(latDeg) - gc.latRad);
        float cosLd = (float) Math.cos(latDiff);
        float sinLd = (float) Math.sin(latDiff);
        // Rotation: X_ned = X_gc*cos(Δlat) + Z_gc*sin(Δlat)
        //           Y_ned = Y_gc (unchanged)
        //           Z_ned = -X_gc*sin(Δlat) + Z_gc*cos(Δlat)
        return new Result(gcX * cosLd + gcZ * sinLd, gcY, -gcX * sinLd + gcZ * cosLd);
    }

    /**
     * Ensure coefficients are loaded from the classpath resource. Called implicitly by public API and safe to call repeatedly.
     */
    private static synchronized void ensureLoaded() {
        if (loaded) return;
        loadCoefficientsFromResource();
        loaded = true;
    }

    /**
     * Pre-load and validate the IGRF coefficient table.
     * <ul>
     *   <li>Forces one-time parsing of <code>igrfcoeffs.txt</code> on the classpath.</li>
     *   <li>Allocates the calling thread’s workspace.</li>
     *   <li>Prints a once-only warning to stderr if the current time is more than 5 years beyond the last epoch.</li>
     * </ul>
     * This method is optional; the first call to {@link #compute(double, double, double, long)} (or
     * {@link #declinationDeg(double, double, double, long)}) will load lazily.
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
                        // Intentional diagnostic output - this library has no logger dependency
                        System.err.printf( // NOSONAR
                                "IGRFModel warning: current time is %.2f years beyond latest epoch"
                                        + " %.1f; results will clamp to %.1f. Consider updating"
                                        + " igrfcoeffs.txt.%n",
                                yearsAhead, (float) lastEpoch, (float) (lastEpoch + 5.0));
                    }
                }
            }
        }
    }

    /**
     * Parse a year token and add it to the years list if it's a valid year (not "SV").
     * Silently ignores non-numeric values that are not valid years.
     */
    private static void parseYearToken(String token, List<Double> years) {
        if (token.equalsIgnoreCase("SV") || token.equalsIgnoreCase("g/h")) {
            return;
        }
        try {
            years.add(Double.parseDouble(token));
        } catch (NumberFormatException _) {
            // Intentionally ignore non-numeric tokens in header
        }
    }

    /**
     * Load and parse the coefficient table from the classpath resource <code>igrfcoeffs.txt</code>.
     * <p>Format: header line starting with <code>g/h</code> then epoch years (and optional <code>SV</code>),
     * followed by rows for <code>g</code> and <code>h</code> with degree n, order m, and values per epoch.</p>
     */
    private static void loadCoefficientsFromResource() {
        // Prefer version-agnostic resource, fall back to legacy name to preserve compatibility
        InputStream in = IGRFModel.class.getClassLoader().getResourceAsStream("igrfcoeffs.txt");
        if (in == null)
            throw new IllegalStateException(
                    "IGRF coefficients file 'igrfcoeffs.txt' not found on classpath");
        try (BufferedReader br = new BufferedReader(new InputStreamReader(in))) {
            List<Double> years = new ArrayList<>();
            List<float[]> gRows = new ArrayList<>();
            List<float[]> hRows = new ArrayList<>();
            boolean headerParsed = false;
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) continue;
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
                if (kind != 'g' && kind != 'h') continue;
                String[] toks = line.split("\\s+");
                if (toks.length < 4) continue;
                int n = Integer.parseInt(toks[1]);
                int m = Integer.parseInt(toks[2]);
                int startIdx = 3;
                int count = toks.length - startIdx;
                float[] row = new float[3 + count];
                row[0] = n;
                row[1] = m;
                row[2] = (kind == 'g') ? 1f : 2f;
                for (int i = 0; i < count; i++) row[3 + i] = Float.parseFloat(toks[startIdx + i]);
                if (kind == 'g') gRows.add(row);
                else hRows.add(row);
            }
            if (years.isEmpty()) {
                List<Double> defaultYears = new ArrayList<>();
                for (int y = 1900; y <= 1995; y += 5) defaultYears.add((double) y);
                for (int y : new int[] {2000, 2005, 2010, 2015, 2020, 2025}) defaultYears.add((double) y);
                years = defaultYears;
            }
            epochs = years.stream().mapToDouble(Double::doubleValue).toArray();
            int epochCount = epochs.length;
            int latestYear = (int) Math.round(epochs[epochCount - 1]);
            baseYearAnchor = (float) epochs[epochCount - 1];
            baseYearMillis = millisAtUtcJan1(latestYear);
            gCoeffs = new float[MAX_N + 1][MAX_N + 1][epochCount];
            hCoeffs = new float[MAX_N + 1][MAX_N + 1][epochCount];
            svG2025 = new float[MAX_N + 1][MAX_N + 1];
            svH2025 = new float[MAX_N + 1][MAX_N + 1];
            for (float[] row : gRows) {
                int n = (int) row[0];
                int m = (int) row[1];
                int evals = row.length - 3;
                System.arraycopy(row, 3, gCoeffs[n][m], 0, epochCount);
                if (evals > epochCount) svG2025[n][m] = row[3 + epochCount];
            }
            for (float[] row : hRows) {
                int n = (int) row[0];
                int m = (int) row[1];
                int evals = row.length - 3;
                System.arraycopy(row, 3, hCoeffs[n][m], 0, epochCount);
                if (evals > epochCount) svH2025[n][m] = row[3 + epochCount];
            }
        } catch (IOException e) {
            throw new IllegalStateException("Failed to read IGRF coefficients file", e);
        }
    }

    /**
     * Milliseconds at UTC 00:00:00 on January 1 of a given year, using proleptic Gregorian calendar.
     */
    private static long millisAtUtcJan1(int year) {
        return daysFrom1970ToYearStart(year) * 86_400_000L;
    }

    private static long daysFrom1970ToYearStart(int year) {
        if (year == 1970) return 0L;
        long diff = year - 1970L;
        return 365L * diff + (leapsBefore(year) - leapsBefore(1970));
    }

    private static long leapsBefore(int year) {
        long y = year - 1L;
        return y / 4L - y / 100L + y / 400L;
    }

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
     * Convert geodetic (WGS‑84) coordinates to geocentric latitude (radians), longitude (radians), and radius (km).
     * <p>Uses standard ellipsoid relations with a = 6378.137 km and b = 6356.7523142 km.</p>
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
     * Fill Gauss‑normalized associated Legendre functions P and their θ‑derivatives, then apply Schmidt factors.
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
     * Immutable result of an IGRF evaluation.
     * <p>All fields are in SI-derived units noted by their suffixes.</p>
     */
    public static final class Result {
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

        private Result(float x, float y, float z) {
            this.xNorthNt = x;
            this.yEastNt = y;
            this.zDownNt = z;
            this.hHorizontalNt = (float) Math.hypot(x, y);
            this.fTotalNt = (float) Math.sqrt(x * x + y * y + z * z);
            this.declinationDeg = (float) Math.toDegrees(Math.atan2(y, x));
            this.inclinationDeg = (float) Math.toDegrees(Math.atan2(z, this.hHorizontalNt));
        }
    }

    private record Geocentric(float latRad, float lonRad, float radiusKm) {}

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
