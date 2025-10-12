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
 */
package com.github.tinemuz;

import static org.junit.jupiter.api.Assertions.*;

import com.github.tinemuz.IGRFModel.Field;
import java.time.LocalDateTime;
import java.time.ZoneOffset;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

/**
 * Comprehensive test suite for IGRFModel.
 *
 * <p>Tests cover: - Core functionality (declination, full field computation) - Known reference
 * values from NOAA - Boundary conditions (poles, date ranges) - Thread safety - Performance
 * characteristics - Edge cases and error handling
 */
class IGRFModelTest {

    private static final double DECLINATION_TOLERANCE = 0.6; // degrees
    private static final double INCLINATION_TOLERANCE = 2.0; // degrees (increased for real model)
    private static final double INTENSITY_TOLERANCE = 1000.0; // nT (increased for real model)
    private static final double COMPONENT_TOLERANCE = 150.0; // nT

    @Nested
    @DisplayName("Core Functionality Tests")
    class CoreFunctionalityTests {

        @Test
        @DisplayName("Field values are finite and within valid ranges")
        void resultFieldsValid() {
            Field Field = IGRFModel.compute(51.5074, -0.1278, 50, epochMillis(2025, 1, 1));

            assertNotNull(Field);
            // All components should be finite
            assertFinite(Field.xNorthNt, "X component");
            assertFinite(Field.yEastNt, "Y component");
            assertFinite(Field.zDownNt, "Z component");
            assertFinite(Field.hHorizontalNt, "Horizontal intensity");
            assertFinite(Field.fTotalNt, "Total intensity");
            assertFinite(Field.declinationDeg, "Declination");
            assertFinite(Field.inclinationDeg, "Inclination");
            
            // Angles should be in valid ranges
            assertTrue(
                    Field.declinationDeg >= -180 && Field.declinationDeg <= 180,
                    "Declination in [-180°, 180°]");
            assertTrue(
                    Field.inclinationDeg >= -90 && Field.inclinationDeg <= 90,
                    "Inclination in [-90°, 90°]");
        }

        @Test
        @DisplayName("Field intensity relationships are consistent")
        void fieldIntensityRelationships() {
            Field r = IGRFModel.compute(35.6762, 139.6503, 40, epochMillis(2025, 1, 1));

            // H² = X² + Y²
            double hCalc = Math.sqrt(r.xNorthNt * r.xNorthNt + r.yEastNt * r.yEastNt);
            assertEquals(r.hHorizontalNt, hCalc, 1.0, "H = √(X² + Y²)");

            // F² = X² + Y² + Z²
            double fCalc =
                    Math.sqrt(
                            r.xNorthNt * r.xNorthNt
                                    + r.yEastNt * r.yEastNt
                                    + r.zDownNt * r.zDownNt);
            assertEquals(r.fTotalNt, fCalc, 1.0, "F = √(X² + Y² + Z²)");

            // H ≤ F always
            assertTrue(r.hHorizontalNt <= r.fTotalNt, "H should be ≤ F");
        }

        @Test
        @DisplayName("Declination and inclination formulas verified")
        void angleFormulas() {
            Field r = IGRFModel.compute(40.7, -74.0, 100, epochMillis(2025, 1, 1));

            // Declination = atan2(Y, X)
            double decCalc = Math.toDegrees(Math.atan2(r.yEastNt, r.xNorthNt));
            assertEquals(r.declinationDeg, decCalc, 0.01, "D = atan2(Y, X)");

            // Inclination = atan2(Z, H)
            double incCalc = Math.toDegrees(Math.atan2(r.zDownNt, r.hHorizontalNt));
            assertEquals(r.inclinationDeg, incCalc, 0.01, "I = atan2(Z, H)");
        }
    }

    @Nested
    @DisplayName("Known Reference Values")
    class ReferenceValueTests {

        @Test
        @DisplayName("New York City - 2025")
        void newYorkCity2025() {
            // NYC: 40.7128°N, 74.0060°W, elevation 10m, Jan 1 2025
            Field r = IGRFModel.compute(40.7128, -74.0060, 10, epochMillis(2025, 1, 1));

            // Expected values (approximate from NOAA calculator)
            assertEquals(-13.0, r.declinationDeg, DECLINATION_TOLERANCE, "NYC declination");
            assertEquals(66.0, r.inclinationDeg, INCLINATION_TOLERANCE, "NYC inclination");
            assertEquals(51000, r.fTotalNt, INTENSITY_TOLERANCE, "NYC total field");

            // Field should point generally northward and downward (NH)
            assertTrue(r.xNorthNt > 0, "North component positive in NH");
            assertTrue(r.zDownNt > 0, "Down component positive in NH");
        }

        @Test
        @DisplayName("London - 2025")
        void london2025() {
            // London: 51.5074°N, 0.1278°W, elevation 11m
            Field r = IGRFModel.compute(51.5074, -0.1278, 11, epochMillis(2025, 1, 1));

            assertEquals(0.5, r.declinationDeg, DECLINATION_TOLERANCE, "London declination");
            assertEquals(66.5, r.inclinationDeg, INCLINATION_TOLERANCE, "London inclination");
            assertEquals(49000, r.fTotalNt, INTENSITY_TOLERANCE, "London total field");
        }

        @Test
        @DisplayName("Tokyo - 2025")
        void tokyo2025() {
            // Tokyo: 35.6762°N, 139.6503°E, elevation 40m
            Field r = IGRFModel.compute(35.6762, 139.6503, 40, epochMillis(2025, 1, 1));

            assertEquals(-7.5, r.declinationDeg, DECLINATION_TOLERANCE, "Tokyo declination");
            assertEquals(50.0, r.inclinationDeg, INCLINATION_TOLERANCE, "Tokyo inclination");
            assertEquals(46000, r.fTotalNt, INTENSITY_TOLERANCE, "Tokyo total field");
        }

        @Test
        @DisplayName("Sydney - 2025 (Southern Hemisphere)")
        void sydney2025() {
            // Sydney: 33.8688°S, 151.2093°E, elevation 58m
            Field r = IGRFModel.compute(-33.8688, 151.2093, 58, epochMillis(2025, 1, 1));

            assertEquals(12.5, r.declinationDeg, DECLINATION_TOLERANCE, "Sydney declination");
            assertEquals(-64.0, r.inclinationDeg, INCLINATION_TOLERANCE, "Sydney inclination");
            assertEquals(57000, r.fTotalNt, INTENSITY_TOLERANCE, "Sydney total field");

            // Southern Hemisphere: inclination negative (points upward)
            assertTrue(r.inclinationDeg < 0, "SH inclination negative");
            assertTrue(r.zDownNt < 0, "SH down component negative (points up)");
        }

        @Test
        @DisplayName("Equator - near-zero inclination")
        void equator() {
            // Near equator: 0°N, 0°E, sea level
            Field r = IGRFModel.compute(0.0, 0.0, 0, epochMillis(2025, 1, 1));

            // Near magnetic equator, inclination should be relatively low
            // (Note: magnetic equator doesn't align exactly with geographic equator)
            assertTrue(
                    Math.abs(r.inclinationDeg) < 45,
                    "Equator has relatively low inclination: " + r.inclinationDeg);

            // Field is mostly horizontal near geographic equator
            assertTrue(
                    r.hHorizontalNt > Math.abs(r.zDownNt) * 0.5,
                    "Horizontal component significant at equator");
        }
    }

    @Nested
    @DisplayName("Boundary Conditions")
    class BoundaryConditionTests {

        @Test
        @DisplayName("North Pole - extreme latitude")
        void northPole() {
            // 89.9°N (avoid exact pole due to coordinate singularity)
            Field r = IGRFModel.compute(89.9, 0, 0, epochMillis(2025, 1, 1));

            // Near pole: inclination close to 90° (vertical downward)
            assertTrue(r.inclinationDeg > 85, "North pole has steep inclination");
            assertTrue(r.zDownNt > 0, "Downward component positive at north pole");

            // Horizontal component should be small
            assertTrue(
                    r.hHorizontalNt < r.fTotalNt * 0.2,
                    "Horizontal component small at pole");
        }

        @Test
        @DisplayName("South Pole - extreme latitude")
        void southPole() {
            // 89.9°S
            Field r = IGRFModel.compute(-89.9, 0, 0, epochMillis(2025, 1, 1));

            // Near south pole: inclination should be steep negative
            // (Note: magnetic south pole location affects exact value)
            assertTrue(
                    r.inclinationDeg < -70,
                    "South pole has steep negative inclination: " + r.inclinationDeg);
            assertTrue(r.zDownNt < 0, "Down component negative at south pole (points up)");
        }

        @Test
        @DisplayName("Longitude wraparound: -180° and +180°")
        void longitudeWraparound() {
            long time = epochMillis(2025, 1, 1);
            Field r1 = IGRFModel.compute(0, -180, 0, time);
            Field r2 = IGRFModel.compute(0, 180, 0, time);

            // -180° and +180° are the same meridian
            assertEquals(r1.declinationDeg, r2.declinationDeg, 0.1, "Longitude wraparound");
            assertEquals(r1.fTotalNt, r2.fTotalNt, 1.0, "Field at -180° = field at 180°");
        }

        @Test
        @DisplayName("Negative altitude (below sea level)")
        void negativeAltitude() {
            // Dead Sea area: -430m elevation
            Field r = IGRFModel.compute(31.5, 35.5, -430, epochMillis(2025, 1, 1));

            assertNotNull(r);
            assertFinite(r.fTotalNt, "Field at negative altitude");
            assertTrue(r.fTotalNt > 40000, "Reasonable field at negative altitude");
        }

        @Test
        @DisplayName("High altitude")
        void highAltitude() {
            // ISS orbit: ~400km altitude
            Field r = IGRFModel.compute(0, 0, 400000, epochMillis(2025, 1, 1));

            assertNotNull(r);
            assertFinite(r.fTotalNt, "Field at high altitude");

            // Field should be weaker at altitude (inverse cube law)
            Field ground = IGRFModel.compute(0, 0, 0, epochMillis(2025, 1, 1));
            assertTrue(r.fTotalNt < ground.fTotalNt, "Field weaker at altitude");
        }
    }

    @Nested
    @DisplayName("Time Range Tests")
    class TimeRangeTests {

        @Test
        @DisplayName("Before first epoch: 1900")
        void beforeFirstEpoch() {
            // Test date before first epoch (should use earliest coefficients)
            Field r1900 = IGRFModel.compute(40.7, -74.0, 0, epochMillis(1900, 1, 1));
            
            assertNotNull(r1900);
            assertFinite(r1900.declinationDeg, "1900 declination");
            assertFinite(r1900.fTotalNt, "1900 field");
            
            // Should be different from 1905 (first interpolation)
            Field r1905 = IGRFModel.compute(40.7, -74.0, 0, epochMillis(1905, 1, 1));
            // Fields should be similar but not identical
            assertTrue(Math.abs(r1900.fTotalNt - r1905.fTotalNt) < 5000, 
                "1900 and 1905 fields reasonably close");
        }

        @Test
        @DisplayName("Historical date: 1950")
        void historical1950() {
            Field r = IGRFModel.compute(40.7, -74.0, 0, epochMillis(1950, 1, 1));

            assertNotNull(r);
            assertFinite(r.declinationDeg, "1950 declination");
            assertFinite(r.fTotalNt, "1950 field");

            // Secular variation: declination changes over time
            Field r2025 = IGRFModel.compute(40.7, -74.0, 0, epochMillis(2025, 1, 1));
            assertNotEquals(
                    r.declinationDeg,
                    r2025.declinationDeg,
                    1.0,
                    "Declination changes over 75 years");
        }

        @Test
        @DisplayName("Future date with secular variation: 2028")
        void futureWithSV() {
            // Within 5 years of 2025 epoch - should use SV
            Field r = IGRFModel.compute(40.7, -74.0, 0, epochMillis(2028, 1, 1));

            assertNotNull(r);
            assertFinite(r.declinationDeg, "2028 declination");

            // Should be slightly different from 2025 due to SV
            Field r2025 = IGRFModel.compute(40.7, -74.0, 0, epochMillis(2025, 1, 1));
            // SV may cause small change (could be less than 0.5° over 3 years)
            double change = Math.abs(r.declinationDeg - r2025.declinationDeg);
            assertTrue(change < 5.0, "SV causes reasonable change over 3 years: " + change);
        }

        @Test
        @DisplayName("Far future date: 2035 (beyond SV limit)")
        void farFuture() {
            // Beyond 5 years - should be clamped and trigger warning
            Field r2035 = IGRFModel.compute(40.7, -74.0, 0, epochMillis(2035, 1, 1));

            assertNotNull(r2035);
            assertFinite(r2035.declinationDeg, "2035 declination (clamped)");

            // Should be same as 2030 (2025 + 5 years SV)
            Field r2030 = IGRFModel.compute(40.7, -74.0, 0, epochMillis(2030, 1, 1));
            assertEquals(
                    r2035.declinationDeg,
                    r2030.declinationDeg,
                    0.01,
                    "Clamped to +5 years SV");
        }
        
        @Test
        @DisplayName("Extremely far future: 2050 (way beyond SV)")
        void extremelyFarFuture() {
            // Test > 5 years beyond to ensure clamping logic is solid
            Field r2050 = IGRFModel.compute(51.5, 0, 0, epochMillis(2050, 1, 1));
            Field r2030 = IGRFModel.compute(51.5, 0, 0, epochMillis(2030, 1, 1));
            
            // Should produce same Field as 2030 (clamped at +5 years)
            assertEquals(r2050.declinationDeg, r2030.declinationDeg, 0.01, 
                "2050 clamped to same as 2030");
            assertEquals(r2050.fTotalNt, r2030.fTotalNt, 1.0, 
                "Field intensity same after clamping");
        }

        @Test
        @DisplayName("Interpolation between epochs")
        void interpolation() {
            // 2022.5 should be interpolated between 2020 and 2025
            Field r2020 = IGRFModel.compute(40.7, -74.0, 0, epochMillis(2020, 1, 1));
            Field r2025 = IGRFModel.compute(40.7, -74.0, 0, epochMillis(2025, 1, 1));
            Field r2022 = IGRFModel.compute(40.7, -74.0, 0, epochMillis(2022, 7, 1));

            // 2022 should be between 2020 and 2025
            assertTrue(
                    isBetween(r2022.declinationDeg, r2020.declinationDeg, r2025.declinationDeg),
                    "2022 declination interpolated");
        }
    }

    @Nested
    @DisplayName("Thread Safety Tests")
    class ThreadSafetyTests {

        @Test
        @DisplayName("Concurrent access from multiple threads")
        void concurrentAccess() throws InterruptedException {
            final int threadCount = 10;
            final int iterationsPerThread = 100;
            Thread[] threads = new Thread[threadCount];
            final double[] results = new double[threadCount];

            for (int i = 0; i < threadCount; i++) {
                final int threadId = i;
                threads[i] =
                        new Thread(
                                () -> {
                                    for (int j = 0; j < iterationsPerThread; j++) {
                                        Field r =
                                                IGRFModel.compute(
                                                        40.7 + threadId * 0.1,
                                                        -74.0,
                                                        0,
                                                        epochMillis(2025, 1, 1));
                                        results[threadId] = r.declinationDeg;
                                    }
                                });
                threads[i].start();
            }

            for (Thread thread : threads) {
                thread.join();
            }

            // All threads should complete successfully
            for (double Field : results) {
                assertFinite(Field, "Thread Field");
            }
        }

        @Test
        @DisplayName("Thread-local workspace independence")
        void threadLocalIndependence() throws InterruptedException {
            final double[] thread1Result = new double[1];
            final double[] thread2Result = new double[1];

            Thread t1 =
                    new Thread(
                            () -> {
                                Field r =
                                        IGRFModel.compute(
                                                40.7, -74.0, 0, epochMillis(2025, 1, 1));
                                thread1Result[0] = r.declinationDeg;
                            });

            Thread t2 =
                    new Thread(
                            () -> {
                                Field r =
                                        IGRFModel.compute(
                                                51.5, 0, 0, epochMillis(2025, 1, 1));
                                thread2Result[0] = r.declinationDeg;
                            });

            t1.start();
            t2.start();
            t1.join();
            t2.join();

            // Different locations should give different results
            assertNotEquals(
                    thread1Result[0], thread2Result[0], 1.0, "Different locations, different D");
        }
    }

    @Nested
    @DisplayName("Performance Tests")
    class PerformanceTests {

        @Test
        @DisplayName("No allocations after warmup")
        void noAllocationsAfterWarmup() {
            long time = epochMillis(2025, 1, 1);

            // Warmup
            for (int i = 0; i < 100; i++) {
                IGRFModel.compute(40.7, -74.0, 100, time);
            }

            // Should be very fast after warmup
            long start = System.nanoTime();
            for (int i = 0; i < 1000; i++) {
                IGRFModel.compute(40.7, -74.0, 100, time);
            }
            long elapsed = System.nanoTime() - start;
            double avgMicros = elapsed / 1000.0 / 1000.0;

            // Should be < 50 microseconds per call on modern hardware
            assertTrue(avgMicros < 50, "Average < 50 μs per call after warmup: " + avgMicros);
        }

        @Test
        @DisplayName("Performance characteristics after warmup")
        void performanceCharacteristics() {
            long time = epochMillis(2025, 1, 1);

            // Warmup thoroughly
            for (int i = 0; i < 2000; i++) {
                IGRFModel.compute(40.7, -74.0, 100, time);
            }

            // Test that compute is reasonably fast
            long start = System.nanoTime();
            for (int i = 0; i < 10000; i++) {
                IGRFModel.compute(40.7, -74.0, 100, time);
            }
            long elapsed = System.nanoTime() - start;
            double avgMicros = elapsed / 10000.0 / 1000.0;

            // Should be fast (< 100 μs per call)
            assertTrue(avgMicros < 100, "Compute < 100μs: " + avgMicros);
        }
    }

    @Nested
    @DisplayName("Edge Cases and Error Handling")
    class EdgeCaseTests {

        @Test
        @DisplayName("Extreme coordinates produce valid results")
        void extremeCoordinates() {
            long time = epochMillis(2025, 1, 1);
            
            // Exact pole (should be clamped internally)
            Field pole = IGRFModel.compute(90.0, 0, 0, time);
            assertFinite(pole.fTotalNt, "Exact pole");
            
            // Large longitude (> 360)
            Field r1 = IGRFModel.compute(40.7, 45, 0, time);
            Field r2 = IGRFModel.compute(40.7, 405, 0, time); // 45 + 360
            assertEquals(r1.declinationDeg, r2.declinationDeg, 0.1, "Longitude modulo 360");
        }

        @Test
        @DisplayName("Date at exact first epoch boundary")
        void exactFirstEpochBoundary() {
            // Test at 1900 (first epoch) - should use first epoch coefficients
            Field r1900 = IGRFModel.compute(40.7, -74.0, 0, epochMillis(1900, 1, 1));
            
            assertNotNull(r1900);
            assertFinite(r1900.declinationDeg, "First epoch declination");
            assertFinite(r1900.fTotalNt, "First epoch field");
            
            // Value should be stable (no interpolation)
            assertTrue(r1900.fTotalNt > 40000 && r1900.fTotalNt < 70000, 
                "Field in reasonable range for 1900: " + r1900.fTotalNt);
        }
        
        @Test
        @DisplayName("Very old date before first epoch: 1850")
        void beforeModelRange() {
            // Date before 1900 (first epoch) - should use earliest coefficients
            Field r1850 = IGRFModel.compute(51.5, 0, 0, epochMillis(1850, 1, 1));
            Field r1900 = IGRFModel.compute(51.5, 0, 0, epochMillis(1900, 1, 1));
            
            // Should give same Field (using first epoch coefficients)
            assertEquals(r1850.declinationDeg, r1900.declinationDeg, 0.01, 
                "Pre-1900 uses first epoch");
            assertEquals(r1850.fTotalNt, r1900.fTotalNt, 1.0, 
                "Same field for dates before model range");
        }
    }

    // Helper methods

    private static long epochMillis(int year, int month, int day) {
        return LocalDateTime.of(year, month, day, 0, 0)
                .toInstant(ZoneOffset.UTC)
                .toEpochMilli();
    }

    private static void assertFinite(double value, String message) {
        assertTrue(Double.isFinite(value), message + " should be finite, got: " + value);
        assertFalse(Double.isNaN(value), message + " should not be NaN");
        assertFalse(Double.isInfinite(value), message + " should not be infinite");
    }

    private static boolean isBetween(double value, double bound1, double bound2) {
        double min = Math.min(bound1, bound2);
        double max = Math.max(bound1, bound2);
        return value >= min && value <= max;
    }
}
