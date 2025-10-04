# IGRF Model Library

A lightweight, high-performance Java library for calculating Earth's magnetic field using the International Geomagnetic Reference Field (IGRF) model. Perfect for navigation applications requiring magnetic declination, inclination, and field intensity calculations.

## Overview

This library provides fast, allocation-free computation of the geomagnetic field vector at any location and time on Earth. It's designed as a minimal dependency suitable for real-time applications, aerospace systems, marine navigation, surveying, and any application requiring conversion between true and magnetic headings.

### Key Features

- **Zero allocation per call** after initialization - ideal for real-time systems
- **Thread-safe** with per-thread workspace caching
- **High accuracy** using the IGRF spherical harmonic model (degree 13)
- **No external dependencies** - pure Java implementation
- **Small footprint** - single class with embedded coefficients
- **WGS-84 ellipsoid** support for geodetic coordinates

## Installation

### Maven

```xml
<dependency>
    <groupId>com.github.tinemuz</groupId>
    <artifactId>igrf-model</artifactId>
    <version>1.0-SNAPSHOT</version>
</dependency>
```

### Gradle

```gradle
implementation 'com.github.tinemuz:igrf-model:1.0-SNAPSHOT'
```

## Quick Start

### Basic Usage - Magnetic Declination

Calculate magnetic declination for converting between true and magnetic headings:

```java
import com.github.tinemuz.IGRFModel;

public class NavigationExample {
    public static void main(String[] args) {
        // Preload coefficients at application startup (recommended)
        IGRFModel.preload();
        
        // Location: New York City (40.7128°N, 74.0060°W), 100m altitude
        double latitude = 40.7128;
        double longitude = -74.0060;
        double altitudeMeters = 100.0;
        long epochMillis = System.currentTimeMillis();
        
        // Get magnetic declination only (fastest method)
        double declination = IGRFModel.declinationDeg(latitude, longitude, altitudeMeters, epochMillis);
        
        System.out.printf("Magnetic declination: %.2f° %s%n", 
            Math.abs(declination), 
            declination < 0 ? "West" : "East");
        
        // Convert true heading to magnetic heading
        double trueHeading = 90.0; // East
        double magneticHeading = trueHeading - declination;
        System.out.printf("True heading %.1f° = Magnetic heading %.1f°%n", 
            trueHeading, magneticHeading);
    }
}
```

### Full Field Computation

Calculate all magnetic field components and derived quantities:

```java
import com.github.tinemuz.IGRFModel;
import com.github.tinemuz.IGRFModel.Result;

public class FullFieldExample {
    public static void main(String[] args) {
        // Preload at startup (one-time operation)
        IGRFModel.preload();
        
        // Location and time
        double latitude = 51.5074;    // London
        double longitude = -0.1278;
        double altitudeMeters = 50.0;
        long epochMillis = System.currentTimeMillis();
        
        // Compute full magnetic field
        Result field = IGRFModel.compute(latitude, longitude, altitudeMeters, epochMillis);
        
        // Access all components
        System.out.printf("North component (X): %.1f nT%n", field.xNorthNt);
        System.out.printf("East component (Y):  %.1f nT%n", field.yEastNt);
        System.out.printf("Down component (Z):  %.1f nT%n", field.zDownNt);
        System.out.printf("Horizontal intensity: %.1f nT%n", field.hHorizontalNt);
        System.out.printf("Total intensity:      %.1f nT%n", field.fTotalNt);
        System.out.printf("Declination:          %.2f°%n", field.declinationDeg);
        System.out.printf("Inclination (dip):    %.2f°%n", field.inclinationDeg);
    }
}
```

## Use Cases

### Aviation Navigation

```java
// Calculate magnetic variation for flight planning
double magneticVariation = IGRFModel.declinationDeg(
    flightPlanLat, flightPlanLon, cruiseAltitude, departureTime);

// Apply to true course to get magnetic course
double magneticCourse = trueCourse - magneticVariation;
```

### Marine Navigation

```java
// Compass deviation correction
Result field = IGRFModel.compute(vesselLat, vesselLon, 0, System.currentTimeMillis());
double magneticDeclination = field.declinationDeg;

// Convert compass heading to true heading
double trueHeading = compassHeading + magneticDeclination + deviationCorrection;
```

### Survey & Geophysics

```java
// Total field intensity for magnetic anomaly detection
Result field = IGRFModel.compute(surveyLat, surveyLon, surveyAlt, surveyTime);
double referenceField = field.fTotalNt;
double measuredField = magnetometer.read();
double anomaly = measuredField - referenceField;
```

### Drone & Autonomous Systems

```java
// Real-time magnetometer calibration reference
public class DroneNavigation {
    static {
        IGRFModel.preload(); // Load once at startup
    }
    
    public void updatePosition(double lat, double lon, double alt) {
        long now = System.currentTimeMillis();
        Result field = IGRFModel.compute(lat, lon, alt, now);
        
        // Use field components for magnetometer calibration
        compassCalibration.setReference(
            field.xNorthNt, field.yEastNt, field.zDownNt);
    }
}
```

## API Reference

### Static Methods

#### `preload()`
```java
public static void preload()
```
Pre-loads IGRF coefficients from the classpath and validates the model. **Recommended to call at application startup** to avoid lazy initialization delays during first use. Prints a warning if the current time is more than 5 years beyond the latest IGRF epoch.

#### `declinationDeg()`
```java
public static double declinationDeg(
    double latDeg, 
    double lonDeg, 
    double altitudeMeters, 
    long epochMillis)
```
Returns only the magnetic declination in degrees. East declination is positive, west is negative. Range: [-180°, +180°).

**Parameters:**
- `latDeg` - Geodetic latitude in degrees (north positive), range typically [-90, +90]
- `lonDeg` - Geodetic longitude in degrees (east positive), any value accepted
- `altitudeMeters` - Altitude above mean sea level in meters (WGS-84)
- `epochMillis` - UTC time in milliseconds since epoch

**Returns:** Declination in degrees east of true north

#### `compute()`
```java
public static Result compute(
    double gdLatitudeDeg, 
    double gdLongitudeDeg, 
    double altitudeMeters, 
    long epochMillis)
```
Computes the full magnetic field vector and derived quantities.

**Returns:** `Result` object containing all field components and values

### Result Class

The `Result` class contains immutable field values:

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `xNorthNt` | float | nT | North component (positive northward) |
| `yEastNt` | float | nT | East component (positive eastward) |
| `zDownNt` | float | nT | Down component (positive downward) |
| `hHorizontalNt` | float | nT | Horizontal intensity √(X² + Y²) |
| `fTotalNt` | float | nT | Total intensity √(X² + Y² + Z²) |
| `declinationDeg` | float | degrees | Declination (variation) = atan2(Y, X) |
| `inclinationDeg` | float | degrees | Inclination (dip) = atan2(Z, H) |

## Mathematical Background

### IGRF Spherical Harmonic Model

The Earth's magnetic field is modeled using spherical harmonics. The magnetic scalar potential $V$ at position $(r, \theta, \phi)$ in spherical coordinates is:

```math
V(r, \theta, \phi) = a \sum_{n=1}^{N} \left(\frac{a}{r}\right)^{n+1} \sum_{m=0}^{n} \left[ g_n^m \cos(m\phi) + h_n^m \sin(m\phi) \right] P_n^m(\cos\theta)
```

Where:
- $a = 6371.2$ km is the reference Earth radius
- $r$ is the geocentric radius
- $\theta$ is the geocentric colatitude (0° at North Pole, 180° at South Pole)
- $\phi$ is the longitude (east positive)
- $g_n^m, h_n^m$ are the Gauss coefficients (Schmidt quasi-normalized)
- $P_n^m$ are the associated Legendre functions (Schmidt quasi-normalized)
- $N = 13$ is the maximum degree and order

### Magnetic Field Components

The magnetic field vector is obtained as the negative gradient of the potential in spherical coordinates:

```math
\mathbf{B} = -\nabla V
```

**Radial component** (positive outward):
```math
B_r = -\frac{\partial V}{\partial r} = \sum_{n=1}^{N} (n+1) \left(\frac{a}{r}\right)^{n+2} \sum_{m=0}^{n} \left[ g_n^m \cos(m\phi) + h_n^m \sin(m\phi) \right] P_n^m(\cos\theta)
```

**Colatitude component** (positive southward):
```math
B_\theta = -\frac{1}{r}\frac{\partial V}{\partial \theta} = -\sum_{n=1}^{N} \left(\frac{a}{r}\right)^{n+2} \sum_{m=0}^{n} \left[ g_n^m \cos(m\phi) + h_n^m \sin(m\phi) \right] \frac{dP_n^m}{d\theta}
```

**Longitude component** (positive eastward):
```math
B_\phi = -\frac{1}{r\sin\theta}\frac{\partial V}{\partial \phi} = \sum_{n=1}^{N} \left(\frac{a}{r}\right)^{n+2} \sum_{m=0}^{n} m \left[ -g_n^m \sin(m\phi) + h_n^m \cos(m\phi) \right] \frac{P_n^m(\cos\theta)}{\sin\theta}
```

### Derived Quantities

**Horizontal Intensity**:
```math
H = \sqrt{X^2 + Y^2}
```

**Total Intensity**:
```math
F = \sqrt{X^2 + Y^2 + Z^2}
```

**Declination** (angle from true north, positive eastward):
```math
D = \arctan\left(\frac{Y}{X}\right)
```

**Inclination** (dip angle, positive downward):
```math
I = \arctan\left(\frac{Z}{H}\right)
```

### Geodetic to Geocentric Conversion

Input coordinates are geodetic (WGS-84). The conversion to geocentric coordinates uses the WGS-84 ellipsoid parameters:
- Semi-major axis: $a = 6378.137$ km
- Semi-minor axis: $b = 6356.7523142$ km
- Flattening: $f = (a-b)/a \approx 1/298.257223563$

**Geocentric latitude** $\phi_c$ from geodetic latitude $\phi_g$:
```math
\phi_c = \arctan\left(\frac{b^2}{a^2} \tan\phi_g\right)
```

**Geocentric radius**:
```math
r = \sqrt{(a^2\cos^2\phi_g + b^2\sin^2\phi_g) + 2h\sqrt{a^2\cos^2\phi_g + b^2\sin^2\phi_g} + h^2}
```

Where $h$ is the altitude above the ellipsoid.

### Time Interpolation

**Between epochs** ($t$ between year $y_i$ and $y_{i+1}$):
```math
g_n^m(t) = g_n^m(y_i) + \frac{t - y_i}{y_{i+1} - y_i}[g_n^m(y_{i+1}) - g_n^m(y_i)]
```

**Beyond last epoch** (using secular variation $\dot{g}_n^m$ for $n \leq 8$):
```math
g_n^m(t) = g_n^m(y_{last}) + \min(\Delta t, 5) \cdot \dot{g}_n^m
```

Where $\Delta t = t - y_{last}$ and is clamped to 5 years maximum.

### Schmidt Quasi-Normalization

The IGRF model uses Schmidt quasi-normalized associated Legendre functions. The Schmidt normalization factor for $P_n^m$ is:

```math
S_n^m = \sqrt{\frac{(n-m)!}{(n+m)!} \cdot (2-\delta_{0m})}
```

Where $\delta_{0m}$ is the Kronecker delta (1 if $m=0$, 0 otherwise).

**Recursive computation** of Schmidt factors:
```math
S_n^0 = S_{n-1}^0 \cdot \frac{2n-1}{n}
```

```math
S_n^m = S_n^{m-1} \cdot \sqrt{\frac{(n-m+1)(2-\delta_{1m})}{n+m}} \quad \text{for } m > 0
```

### Associated Legendre Functions

The library computes **Gauss-normalized** associated Legendre functions using stable recurrence relations, then multiplies by Schmidt factors.

**Initial values**:
```math
P_0^0(\cos\theta) = 1
```

**Diagonal recurrence** ($n = m$):
```math
P_n^n(\cos\theta) = \sin\theta \cdot P_{n-1}^{n-1}(\cos\theta)
```

**Vertical recurrence** ($m = n-1$):
```math
P_n^{n-1}(\cos\theta) = \cos\theta \cdot P_{n-1}^{n-1}(\cos\theta)
```

**General recurrence** ($m < n-1$):
```math
P_n^m(\cos\theta) = \cos\theta \cdot P_{n-1}^m(\cos\theta) - k_{nm} \cdot P_{n-2}^m(\cos\theta)
```

Where:
```math
k_{nm} = \frac{(n-1)^2 - m^2}{(2n-1)(2n-3)}
```

**Derivative recurrence** for $\frac{dP_n^m}{d\theta}$:
```math
\frac{dP_n^m}{d\theta} = -\sin\theta \cdot P_{n-1}^m + \cos\theta \cdot \frac{dP_{n-1}^m}{d\theta} - k_{nm} \cdot \frac{dP_{n-2}^m}{d\theta}
```

### Computational Optimizations

The library uses several optimizations to achieve zero-allocation, high-performance computation:

**1. Trigonometric Series Caching**

Instead of computing $\sin(m\phi)$ and $\cos(m\phi)$ for each order $m$, we use angle addition formulas:

```math
\sin((m+1)\phi) = \sin(m\phi)\cos(\phi) + \cos(m\phi)\sin(\phi)
```

```math
\cos((m+1)\phi) = \cos(m\phi)\cos(\phi) - \sin(m\phi)\sin(\phi)
```

Starting with $\sin(\phi)$ and $\cos(\phi)$, all higher orders are computed via recurrence.

**2. Radius Power Caching**

Powers of $(a/r)$ are cached:

```math
\left(\frac{a}{r}\right)^{n+1} = \left(\frac{a}{r}\right)^n \cdot \left(\frac{a}{r}\right)
```

This avoids repeated exponentiation operations.

**3. Workspace Caching**

A thread-local workspace stores:
- Associated Legendre functions $P_n^m(\cos\theta)$ for current $\theta$
- Derivatives $\frac{dP_n^m}{d\theta}$
- Schmidt-normalized versions
- Trigonometric series for current $\phi$
- Radius powers for current $r$

Values are only recomputed when coordinates change.

### Geocentric to Geodetic Rotation

After computing field components in the geocentric frame $(B_r, B_\theta, B_\phi)$, we rotate to the local geodetic NED frame. The rotation angle is the difference between geodetic and geocentric latitudes:

```math
\Delta\phi = \phi_{\text{geodetic}} - \phi_{\text{geocentric}}
```

**Rotation to NED frame**:
```math
\begin{bmatrix} X_{\text{north}} \\ Y_{\text{east}} \\ Z_{\text{down}} \end{bmatrix} = 
\begin{bmatrix} \cos(\Delta\phi) & 0 & \sin(\Delta\phi) \\ 0 & 1 & 0 \\ -\sin(\Delta\phi) & 0 & \cos(\Delta\phi) \end{bmatrix}
\begin{bmatrix} B_\theta \\ B_\phi \\ -B_r \end{bmatrix}
```

Note: $B_\theta$ points south (increasing $\theta$), so it becomes the north component. $B_r$ points outward, so its negative becomes the down component.

### Heading Conversions

**True to Magnetic**:
```math
\text{Heading}_{\text{mag}} = \text{Heading}_{\text{true}} - D
```

**Magnetic to True**:
```math
\text{Heading}_{\text{true}} = \text{Heading}_{\text{mag}} + D
```

Where $D$ is the declination (variation).

## Coordinate Systems

### Input Coordinates (Geodetic - WGS-84)
- **Latitude**: -90° (South Pole) to +90° (North Pole)
- **Longitude**: -180° to +180° or 0° to 360° (both accepted)
- **Altitude**: Meters above WGS-84 ellipsoid (negative values allowed)

### Output Frame (NED - North-East-Down)
- **X**: Points to geographic north
- **Y**: Points to geographic east  
- **Z**: Points downward toward Earth's center

All output components are in nanoteslas (nT).

### Typical Field Values

The Earth's magnetic field varies by location:

| Location | Total Field (nT) | Declination | Inclination |
|----------|------------------|-------------|-------------|
| **Equator** | ~25,000-40,000 nT | ±20° | ~0° (horizontal) |
| **Mid-latitudes** | ~45,000-55,000 nT | ±30° | 50-70° |
| **Polar regions** | ~55,000-65,000 nT | Variable | ~90° (vertical) |

**Example values (2025)**:
- New York (40.7°N, 74.0°W): F ≈ 52,000 nT, D ≈ -13° W, I ≈ 67°
- London (51.5°N, 0.1°W): F ≈ 49,000 nT, D ≈ 0° to +1° E, I ≈ 66°
- Sydney (33.9°S, 151.2°E): F ≈ 57,000 nT, D ≈ +12° E, I ≈ -64°
- Tokyo (35.7°N, 139.7°E): F ≈ 46,000 nT, D ≈ -7° W, I ≈ 50°

### Model Accuracy

The IGRF-13 model provides the following typical accuracies for 2015-2025:

| Quantity | Accuracy |
|----------|----------|
| **Declination (D)** | ±0.5° globally, ±1° near magnetic poles |
| **Inclination (I)** | ±0.5° globally, ±1° near magnetic equator |
| **Total Intensity (F)** | ±100-200 nT globally |
| **Horizontal Components (X, Y)** | ±100-150 nT |
| **Vertical Component (Z)** | ±150-250 nT |

**Important notes**:
- Accuracy degrades beyond 5 years past the latest epoch
- Local magnetic anomalies can cause deviations of thousands of nT
- For precise survey work, consider regional models or local magnetic surveys
- Update coefficients regularly for best results

## Time Handling

The library uses epoch milliseconds (UTC) and internally converts to decimal years:

- **Historical dates**: Linear interpolation between IGRF epochs
- **Recent dates**: Uses latest epoch coefficients
- **Future dates**: Applies secular variation (SV) up to +5 years beyond last epoch
- **Far future**: Clamps to +5 years with a one-time warning

### Getting Epoch Milliseconds

```java
// Current time
long now = System.currentTimeMillis();

// Specific date/time
LocalDateTime dt = LocalDateTime.of(2025, 1, 15, 12, 0);
long epochMillis = dt.toInstant(ZoneOffset.UTC).toEpochMilli();

// From java.util.Date
long epochMillis = date.getTime();
```

## Performance Characteristics

### Initialization
- **First call**: ~50-100ms (loads coefficients from classpath)
- **Subsequent calls**: No additional loading
- **Memory**: ~200KB for coefficients + thread-local workspace (~20KB per thread)

### Runtime Performance
- **Single evaluation**: ~5-15 microseconds (after warmup)
- **Allocations**: Zero per call after first use
- **Thread safety**: Lock-free via ThreadLocal workspace

### Benchmark Example
```java
// Warmup
for (int i = 0; i < 1000; i++) {
    IGRFModel.declinationDeg(40.7, -74.0, 100, System.currentTimeMillis());
}

// Benchmark
long start = System.nanoTime();
for (int i = 0; i < 100_000; i++) {
    IGRFModel.declinationDeg(40.7, -74.0, 100, System.currentTimeMillis());
}
long elapsed = System.nanoTime() - start;
System.out.printf("Average: %.2f μs per call%n", elapsed / 100_000.0 / 1000.0);
```

## Technical Details

### IGRF Model
- **Spherical harmonic expansion** to degree and order 13
- **Schmidt quasi-normalized** associated Legendre functions
- **Reference radius**: 6371.2 km
- **Epoch coverage**: 1900-2025+ with secular variation

### Computational Approach
1. Convert WGS-84 geodetic coordinates to geocentric (latitude, radius)
2. Compute associated Legendre functions via stable recurrence relations
3. Apply Schmidt normalization factors
4. Evaluate spherical harmonic expansion
5. Compute field gradients in geocentric frame
6. Rotate back to geodetic NED frame

### Caching Strategy
- **Per-thread workspace**: Avoids allocation and contention
- **Cached values**: Legendre functions, trigonometric series, radius powers
- **Invalidation**: Automatic when coordinates change

## Common Questions

### What is magnetic declination?
Magnetic declination (or variation) is the angle between true north (geographic north) and magnetic north. It varies by location and changes over time due to the Earth's dynamic magnetic field.

### How do I convert between true and magnetic headings?

```java
// True to Magnetic: subtract declination
double magneticHeading = trueHeading - declination;

// Magnetic to True: add declination
double trueHeading = magneticHeading + declination;
```

**Example**: At a location with 10° East declination:
- True heading 90° (East) = Magnetic heading 80°
- Magnetic heading 270° (West) = True heading 280°

### What about compass deviation?
This library provides magnetic declination (variation), which is the difference between true and magnetic north. Compass deviation (caused by local magnetic interference) must be determined separately through compass calibration and applied in addition to declination.

### How accurate is the model?
The IGRF model provides accuracy of:
- **Declination**: ~0.5° or better for most locations
- **Total intensity**: ~100-200 nT globally

Accuracy degrades beyond 5 years past the latest epoch. Update coefficients regularly for best results.

### Can I use this for historical data?
Yes, the model covers 1900-2025 with linear interpolation between epochs. For precise historical work, ensure you have coefficients for the time period of interest.

### Thread safety?
Fully thread-safe. Each thread gets its own workspace cache automatically via ThreadLocal storage.

## Updating IGRF Coefficients

The library reads coefficients from `igrfcoeffs.txt` on the classpath. To update:

1. Download latest IGRF coefficients from [NOAA](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)
2. Replace `src/main/resources/igrfcoeffs.txt`
3. Rebuild and redeploy

The file format is the standard IGRF table with header row starting with `g/h` followed by epoch years and `SV` column.

## Example Applications

### Flight Management System
```java
public class FlightComputer {
    private final IGRFModel model;
    
    public FlightComputer() {
        IGRFModel.preload();
    }
    
    public NavigationData getNavigationData(Position pos, long time) {
        double declination = IGRFModel.declinationDeg(
            pos.latitude, pos.longitude, pos.altitude, time);
        
        return new NavigationData()
            .withMagneticVariation(declination)
            .withTrueTrack(calculateTrueTrack())
            .withMagneticTrack(calculateTrueTrack() - declination);
    }
}
```

### Magnetometer Calibration
```java
public class MagnetometerCalibrator {
    public void calibrate(SensorData sensor, Position pos) {
        Result field = IGRFModel.compute(
            pos.lat, pos.lon, pos.alt, System.currentTimeMillis());
        
        // Expected field vector in NED frame
        Vector3D expected = new Vector3D(
            field.xNorthNt, field.yEastNt, field.zDownNt);
        
        // Measured field from magnetometer
        Vector3D measured = sensor.getMagneticField();
        
        // Compute calibration parameters
        CalibrationMatrix cal = computeCalibration(expected, measured);
        sensor.applyCalibration(cal);
    }
}
```

## Troubleshooting

### Warning: "requested time is X years beyond latest epoch"
Your computation time is more than 5 years past the latest IGRF epoch in the coefficients file. Results are clamped to +5 years with secular variation. Update your `igrfcoeffs.txt` file with the latest IGRF model.

### IllegalStateException: "IGRF coefficients file not found"
The `igrfcoeffs.txt` file is missing from the classpath. Ensure it's in `src/main/resources/` and properly packaged in your JAR.

### Performance degradation
If you notice performance issues:
1. Call `preload()` at startup to avoid lazy initialization
2. Reuse the same thread for computations when possible (to benefit from cached workspace)
3. Avoid creating excessive short-lived threads

## Dependencies

**None** - This is a standalone library with no external dependencies.

## License

[Specify your license here]

## References

- [IGRF Model Documentation](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)
- [NOAA Geomagnetic Calculators](https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml)
- [WGS-84 Specification](https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84)

## Contributing

Contributions are welcome! Please ensure:
- Code follows existing style conventions
- All tests pass
- Documentation is updated for API changes
- SonarQube quality gates are maintained

## Version History

### 1.0-SNAPSHOT (Current)
- Initial release
- IGRF-14 model support (epochs 1900-2030)
- Zero-allocation design
- Thread-safe implementation
- WGS-84 geodetic coordinate support
